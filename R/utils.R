
#' Make a named bed into a list of GRanges objects to be used with pairwise
#' functions, mostly.
#'
#' @param bedpath A path to a bedfile
#' @importFrom rtracklayer import
#' @return A list of GRanges objects
#' @export
#'
import_named_bed_into_list <- function(bedpath) {
  gr <- rtracklayer::import(bedpath, format = "BED")
  # this is more efficient than unique(name), for obvious reasons
  labels <- levels(factor(gr$name))

  gr_list <- lapply(labels, .filter_granges_by_name, gr = gr)
  names(gr_list) <- labels
  gr_list
}


## Helpers --------

.filter_granges_by_name <- function(gr, name) {
  gr[gr$name == name, ]
}
