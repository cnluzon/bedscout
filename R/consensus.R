
# NOTE: Intentional left out importFrom GenomicRanges reduce, as it writes
# two importFrom on the NAMESPACE that are for reduce functions and I get a
# warning on the documentation check.

#' Find a consensus set of loci across a set of GRanges objects.
#'
#' Ignores strand.
#'
#' @param grlist List of GRanges for loci groups to be compared.
#' @param min_consensus If a locus appears in at least min_consensus loci sets,
#'   it will be kept. min_consensus must be a number between 1 and length(grlist)
#' @param resize Resize the GRanges objects to a fixed size before the check.
#'   If NULL or 0, the loci are not resized.
#' @param anchor If resize, where to anchor. By default is set to center.
#'
#' @return A GRanges object with the consensus list.
#' @importFrom dplyr filter
#' @importFrom GenomicRanges coverage makeGRangesFromDataFrame resize
#' @importFrom methods as
#' @export
#'
#' @examples
#' gr1 <- GenomicRanges::GRanges(
#'   seqnames = c("chr1"), IRanges::IRanges(c(11, 20), c(15, 30)), strand = "-"
#' )
#' gr2 <- GenomicRanges::GRanges(
#'   seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+"
#' )
#'
#' loci_consensus(list(gr1, gr2), min_consensus = 2)
loci_consensus <- function(grlist, min_consensus = 1, resize = NULL, anchor = "center") {

  .validate_min_consensus(length(grlist), min_consensus)
  if (!is.null(resize)) {
    if (resize > 0) {
      grlist <- lapply(grlist, GenomicRanges::resize, fix = anchor, width = resize)
    }
  }
  cov_gr <- methods::as(GenomicRanges::coverage(do.call(c, unname(grlist))), "GRanges") %>%
    data.frame() %>%
    dplyr::filter(.data$score >= min_consensus) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  # There might be adjacent ranges with different coverage that need to
  # be merged after
  cov_gr <- GenomicRanges::reduce(cov_gr)
  cov_gr
}

## Helpers --------------

.validate_min_consensus <- function(list_length, min_consensus) {
  if (min_consensus < 1) {
    stop(paste("min_consensus must be a positive number:", min_consensus))
  }
  if (min_consensus > list_length) {
    msg <- paste0(
      "min_consensus must be smaller or equal to size of the list (",
      list_length,
      "), got ",
      min_consensus
    )
    stop(msg)
  }
}
