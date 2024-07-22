
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


#' Peaks universe for a set of conditions and replicates
#'
#' @param peak_files A list where each element is a list of peak files
#' @param min_consensus Minimum consensus needed to merge peaks from the same condition (replicates)
#' @param resize_pre Resize replicate peaks before calling consensus - Default NULL
#' @param resize_across Resize peaks before merging across conditions.
#'   Resizing to a small size here will make sure that proximal peaks will not
#'   get merged, but might incur in certain level of overlap if resize_post is
#'   larger.
#' @param resize_post Resize final group of peaks to a fixed size. Default NULL
#' @param anchor Where to anchor for resizing. Defaults to center.
#'
#' @return A GRanges objects with the peaks universe
#' @export
#'
peaks_universe <- function(peak_files, min_consensus, resize_pre, resize_across, resize_post = NULL, anchor = "center") {
  peaks <- lapply(peak_files, .import_peaklist)
  gr_universe(peaks, min_consensus, resize_pre, resize_across, resize_post, anchor)
}


#' GRanges universe for a set of conditions and replicates
#'
#' @param grlist A list where each element is a list of GRanges objects
#' @inheritParams peaks_universe
#' @export
gr_universe <- function(grlist, min_consensus, resize_pre = NULL, resize_across = NULL, resize_post = NULL, anchor = "center") {
  peaks_merged <- lapply(
    grlist,
    loci_consensus,
    min_consensus = min_consensus,
    resize = resize_pre,
    anchor = anchor
  )

  peaks_universe <- loci_consensus(
    peaks_merged,
    min_consensus = 1,
    resize = resize_across,
    anchor = anchor
  )

  .resize_gr(peaks_universe, resize_post, anchor)
}

## Helpers --------------
.import_peaklist <- function(peakfiles) {
  lapply(peakfiles, rtracklayer::import)
}

.resize_gr <- function(gr, size, anchor) {
  result <- gr
  if (!is.null(size)) {
    if (size > 0) {
      result <- GenomicRanges::resize(gr, size, fix = anchor)
    }
  }
  result
}

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
