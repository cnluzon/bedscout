
#' Total overlapping base pairs between two GRanges objects.
#'
#' @param gr1 GRanges object
#' @param gr2 GRanges object
#' @param ignore.strand Whether to ignore strand when calculating overlap.
#'
#' @return Integer representing the total bp overlap between two sets of loci.
#' @importFrom GenomicRanges intersect width
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom IRanges IRanges
#' @export
#'
#' @examples
#' gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
#' gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
#' total_bp_overlap(gr_1, gr_2, ignore.strand = TRUE)
#' total_bp_overlap(gr_1, gr_2, ignore.strand = FALSE)
total_bp_overlap <- function(gr1, gr2, ignore.strand = TRUE) {
  sum(
    width(
      GenomicRanges::intersect(gr1, gr2, ignore.strand = ignore.strand)
    )
  )
}
