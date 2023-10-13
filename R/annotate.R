
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

  # One would think there is a function to expand from the sides of intervals
  expanded_gr <- purrr::reduce(
    list(
      GenomicRanges::flank(gr, distance_cutoff, start = TRUE),
      gr,
      GenomicRanges::flank(gr, distance_cutoff, start = FALSE)
    ),
    GenomicRanges::union,
    ignore.strand = ignore.strand
  )

  nearest_hits <- data.frame(
    findOverlaps(expanded_gr, feat_gr, ignore.strand = ignore.strand)
  )

  feat_df <- data.frame(feat_gr) %>%
    tibble::rowid_to_column("gene_row") %>%
    dplyr::select(dplyr::all_of(c("gene_row", name_field)))

  id_cols <- c("seqnames", "start", "end", "width", "strand")
  cols_select <- c(id_cols, name_field)

  gr_found <- data.frame(gr[nearest_hits$queryHits, ]) %>%
    dplyr::mutate(nearest_id = nearest_hits$subjectHits) %>%
    dplyr::left_join(feat_df, by=c("nearest_id"="gene_row")) %>%
    dplyr::select(dplyr::all_of(cols_select)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(id_cols))) %>%
    dplyr::summarise("nearby_features" = paste(.data[[name_field]], collapse=","))

  makeGRangesFromDataFrame(
    data.frame(gr) %>% # if nearby_name exists in the input gr, drop it
      dplyr::select(-dplyr::any_of("nearby_features")) %>%
      dplyr::left_join(gr_found, by = id_cols),
    keep.extra.columns = TRUE
  )
}

