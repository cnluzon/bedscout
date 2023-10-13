
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
