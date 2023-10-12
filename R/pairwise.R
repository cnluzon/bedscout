
#' Calculate pairwise enrichment matrix between two lists of GRanges objects
#'
#' @param grlist_1 A named list of GRanges
#' @param grlist_2 A named list of GRanges
#' @param genome_size Genome size
#' @param ignore.strand Ignore strand to calculate overlaps.
#'
#' @importFrom purrr map2
#' @return A data frame in long format with the enrichment values calculated
#'   per pair.
#' @export
#'
#' @examples
#' gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
#' gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")
#'
#' gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")
#' gr_b2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 25), strand = "+")
#'
#' gsize <- 30
#'
#' enrichment <- pairwise_enrichment(
#'   list("A" = gr_a1, "B" = gr_a2),
#'   list("X" = gr_b1, "Y" = gr_b2),
#'   genome_size = gsize,
#'   ignore.strand = FALSE
#' )
pairwise_enrichment <- function(grlist_1, grlist_2, genome_size, ignore.strand = TRUE) {
  if (is.null(names(grlist_1)) || is.null(names(grlist_2))) {
    stop("GRanges lists must be named")
  }
  combinations <- expand.grid(names(grlist_1), names(grlist_2))
  names(combinations) <- c("gr1", "gr2")

  result <- mapply(
    fc_enrichment,
    gr1 = grlist_1[combinations$gr1],
    gr2 = grlist_2[combinations$gr2],
    genome_size = genome_size,
    ignore.strand = ignore.strand,
    USE.NAMES = FALSE,
    SIMPLIFY = TRUE
  )

  cbind(combinations, enrichment = result)
}


#' Calculate enrichment combinations of different elements within a list
#'
#' @param grlist Named list of GRanges objects
#' @param genome_size Genome size
#' @param ignore.strand If true, strand is not considered for overlaps.
#'
#' @importFrom utils combn
#' @return A data frame in long format with the enrichment values calculated
#'   per pair.
#' @export
#'
#' @examples
#' gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
#' gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")
#' gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")
#'
#' enrichment <- combinations_enrichment(
#'   list("A" = gr_a1, "B" = gr_a2, "X" = gr_b1),
#'   genome_size = 30,
#'   ignore.strand = FALSE
#' )
combinations_enrichment <- function(grlist, genome_size, ignore.strand = TRUE) {
  if (is.null(names(grlist))) {
    stop("GRanges list must be named")
  }

  combinations <- as.data.frame(t(combn(names(grlist),2, simplify = TRUE)))
  names(combinations) <- c("gr1", "gr2")

  result <- mapply(
    fc_enrichment,
    gr1 = grlist[combinations$gr1],
    gr2 = grlist[combinations$gr2],
    genome_size = genome_size,
    ignore.strand = ignore.strand,
    USE.NAMES = FALSE,
    SIMPLIFY = TRUE
  )

  cbind(combinations, enrichment = result)
}
