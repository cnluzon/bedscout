
#' Total overlapping base pairs between two GRanges objects.
#'
#' @param gr1 GRanges object
#' @param gr2 GRanges object
#' @param ignore.strand Whether to ignore strand when calculating overlap.
#'
#' @return Integer representing the total bp overlap between two sets of loci.
#' @importFrom GenomicRanges intersect width
#' @export
#'
#' @examples
total_bp_overlap <- function(gr1, gr2, ignore.strand = TRUE) {
  sum(
    width(
      GenomicRanges::intersect(gr1, gr2, ignore.strand = ignore.strand)
    )
  )
}
