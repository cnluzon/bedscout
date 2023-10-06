#' Find the intersection of a list of GRanges objects.
#'
#' Note that if only one element is in the list, it returns the same list.
#'
#' @param grlist A list of GRanges
#' @param ignore.strand If TRUE, overlaps are counted without taking strand into
#'   account.
#'
#' @return A GRanges object containing overlapping loci
#' @importFrom purrr reduce
#' @importFrom GenomicRanges intersect
#' @export
#'
#' @examples
#' gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
#' gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
#' gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 18), strand = "+")
#' intersection <- grlist_intersect(list(gr_1, gr_2, gr_3), ignore.strand = TRUE)
grlist_intersect <- function(grlist, ignore.strand = TRUE) {
  purrr::reduce(grlist, GenomicRanges::intersect, ignore.strand = ignore.strand)
}
