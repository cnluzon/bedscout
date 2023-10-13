
#' Annotate each locus in a GRanges with the single best overlap in another
#' feature GRanges. Best overlap is calculated as jaccard index:
#'
#' length(intersect) / length(union)
#'
#' @param gr GRanges to annotate
#' @param feature_gr GRanges with features - Needs to have a name field
#' @param name_field Name of the column that contains the feature names in feature_gr
#' @param minoverlap Minimum overlap to consider.
#' @param ignore.strand Whether to ignore strand on annotation.
#'
#' @importFrom S4Vectors Pairs queryHits
#' @importFrom GenomicRanges findOverlaps pintersect punion makeGRangesFromDataFrame
#' @importFrom IRanges findOverlaps
#' @importFrom dplyr `%>%` .data
#'
#' @return GRanges gr annotated with a "feature" metadata column with the name
#'  of the feature that best overlaps with it, NA if none do.
#' @export
#'
#' @examples
#'
#' features_gr <- GenomicRanges::GRanges(
#'   seqnames = c("chr1", "chr1"),
#'   IRanges::IRanges(c(10,22), c(20,30)),
#'   strand = c("-", "-"),
#'   name = c("Feat_A", "Feat_B")
#' )
#'
#' gr <- GenomicRanges::GRanges(
#'   seqnames = c("chr1"),
#'   IRanges::IRanges(15, 25),
#'   strand = "+"
#' )
#'
#' impute_feature(gr, features_gr, "name", ignore.strand = TRUE)
impute_feature <- function(gr, feature_gr, name_field, minoverlap = 1L, ignore.strand = TRUE) {
  # Again this problem with levels not exactly matching, which happens a lot
  hits <- suppressWarnings(
    findOverlaps(
      gr,
      feature_gr,
      ignore.strand = ignore.strand,
      minoverlap = minoverlap
    )
  )

  p <- Pairs(gr, feature_gr, hits = hits)

  # Again this problem with levels not exactly matching, which happens a lot
  scores <- suppressWarnings(
    width(pintersect(p, ignore.strand = ignore.strand))/width(punion(p, ignore.strand = ignore.strand))
  )

  annotated_hits <- cbind(
      data.frame(hits), scores
    ) %>%
    dplyr::rename(score = 3) %>%
    dplyr::group_by(queryHits) %>%
    dplyr::slice_max(.data$score, na_rm = TRUE)

  best_annotated_df <- cbind(
    cbind(
      data.frame(gr[annotated_hits$queryHits, ]),
      feature = data.frame(feature_gr)[annotated_hits$subjectHits, name_field]
    ), score = annotated_hits$score
  )

  id_cols <- colnames(data.frame(gr))
  makeGRangesFromDataFrame(
    data.frame(gr) %>%
      dplyr::left_join(data.frame(best_annotated_df), by = id_cols),
    keep.extra.columns = TRUE
  )
}


#' Annotates a GRanges object with nearby features in another GRanges object
#'
#' Distance is measured as distance from end and start depending whether the
#' target locus is downstream or upstream.
#'
#' @param gr GRanges object to be annotated
#' @param feat_gr GRanges features to be annotated.
#' @param name_field Name of the column that contains the annotation name in feat_gr
#' @param distance_cutoff Maximum distance around the promoter of the feature
#' @param ignore.strand Ignore strand when annotating. Default FALSE
#'
#' @importFrom tibble rowid_to_column
#' @importFrom IRanges findOverlaps
#' @importFrom dplyr .data `%>%`
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @return A GRanges object with a nearby_features field in it containing names
#'   of overlapping features in that distance range (separated by commas)
#' @export
#'
#' @examples
#'
#' features_gr <- GenomicRanges::GRanges(
#'   seqnames = c("chr1", "chr1"),
#'   IRanges::IRanges(c(10,22), c(20,30)),
#'   strand = c("-", "-"),
#'   name = c("Feat_A", "Feat_B")
#' )
#'
#' gr <- GenomicRanges::GRanges(
#'   seqnames = c("chr1"),
#'   IRanges::IRanges(15, 25),
#'   strand = "+"
#' )
#'
#' annotate_nearby_features(gr, features_gr, "name", ignore.strand = TRUE)
annotate_nearby_features <- function(gr, feat_gr, name_field, distance_cutoff = 2500, ignore.strand = TRUE) {

  # We only want to get overlaps so no need to merge overlapping things, which
  # complicates calculation. Later only unique names are kept
  expanded_feat_gr <- c(
      GenomicRanges::flank(feat_gr, distance_cutoff, start = TRUE),
      feat_gr,
      GenomicRanges::flank(feat_gr, distance_cutoff, start = FALSE)
    )

  annotate_overlapping_features(
    gr,
    expanded_feat_gr,
    name_field,
    ignore.strand = ignore.strand
  )
}


#' Annotates a GRanges object with overlapping features in another GRanges object
#'
#' Distance is measured as distance from end and start depending whether the
#' target locus is downstream or upstream.
#'
#' @param gr GRanges object to be annotated
#' @param feat_gr GRanges features to be annotated.
#' @param name_field Name of the column that contains the annotation name in feat_gr
#' @param minoverlap Minimum overlapping (used by GenomicRanges findOverlaps)
#' @param ignore.strand Ignore strand when annotating. Default FALSE
#'
#' @importFrom tibble rowid_to_column
#' @importFrom IRanges findOverlaps
#' @importFrom dplyr .data `%>%`
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @return A GRanges object with a nearby_features field in it containing names
#'   of overlapping features in that distance range (separated by commas)
#' @export
#'
#' @examples
#'
#' features_gr <- GenomicRanges::GRanges(
#'   seqnames = c("chr1", "chr1"),
#'   IRanges::IRanges(c(10,22), c(20,30)),
#'   strand = c("-", "-"),
#'   name = c("Feat_A", "Feat_B")
#' )
#'
#' gr <- GenomicRanges::GRanges(
#'   seqnames = c("chr1"),
#'   IRanges::IRanges(15, 25),
#'   strand = "+"
#' )
#'
#' annotate_overlapping_features(gr, features_gr, "name", ignore.strand = TRUE)
annotate_overlapping_features <- function(gr, feat_gr, name_field, minoverlap = 1L, ignore.strand = TRUE) {
  # Again this problem with levels not exactly matching, which happens a lot
  hits <- suppressWarnings(
    findOverlaps(
      gr,
      feat_gr,
      ignore.strand = ignore.strand,
      minoverlap = minoverlap
    )
  )

  id_cols <- colnames(data.frame(gr))

  # Take the ranges from query and the names from subject hits, and then
  # comma-collapse the names
  annotated_hits <- cbind(
    data.frame(
      gr[S4Vectors::queryHits(hits), ]),
    name = S4Vectors::mcols(feat_gr[S4Vectors::subjectHits(hits), ])[[name_field]]
  ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) %>%
    dplyr::summarise("nearby_features" = paste(unique(.data[["name"]]), collapse=","))

  # Left-join so the non-hits are kept
  makeGRangesFromDataFrame(
    data.frame(gr) %>%
      dplyr::left_join(data.frame(annotated_hits), by = id_cols),
    keep.extra.columns = TRUE
  )
}
