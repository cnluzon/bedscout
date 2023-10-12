
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
  score_func <- purrr::partial(
    fc_enrichment,
    genome_size = genome_size,
    ignore.strand = ignore.strand
  )
  df <- pairwise_score(grlist_1, grlist_2, score_func)
  names(df)[names(df) == "score"] <- "enrichment"
  df
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
  score_func <- purrr::partial(
    fc_enrichment,
    genome_size = genome_size,
    ignore.strand = ignore.strand
  )
  df <- combinations_score(grlist, score_func)
  names(df)[names(df) == "score"] <- "enrichment"
  df
}


#' Apply a given score function to all pairs in a list of GRanges objects
#'
#' @param grlist named List of GRanges objects
#' @param score_func Any scoring function (must take a pair of GRanges objects)
#'
#' @return A data frame
#' @export
#'
#' @examples
#' gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
#' gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")
#' gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")
#'
#' # Fix the function parameters
#' score_func <- purrr::partial(
#'   fc_enrichment,
#'   genome_size = 30,
#'   ignore.strand = TRUE
#' )
#'
#' combinations_score(
#'   list("A" = gr_a1, "B" = gr_a2, "X" = gr_b1),
#'   score_func
#' )
combinations_score <- function(grlist, score_func) {
  if (is.null(names(grlist))) {
    stop("GRanges list must be named")
  }
  combinations <- as.data.frame(t(combn(names(grlist),2, simplify = TRUE)))
  names(combinations) <- c("gr1", "gr2")
  result <- mapply(
    score_func,
    grlist[combinations$gr1],
    grlist[combinations$gr2],
    USE.NAMES = FALSE,
    SIMPLIFY = TRUE
  )
  cbind(combinations, score = result)
}

#' Apply a given score function to all pairwise combinations of two named
#' GRanges lists
#'
#' @param grlist_1 A named List of GRanges objects
#' @param grlist_2 A named List of GRanges objects
#' @param score_func Any scoring function (must take a pair of GRanges objects)
#'
#' @return A data frame
#' @export
#'
#' @examples
#' gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
#' gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")
#' gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")
#' gr_b2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 25), strand = "+")
#'
#' # Fix the function parameters
#' score_func <- purrr::partial(
#'   fc_enrichment,
#'   genome_size = 30,
#'   ignore.strand = FALSE
#' )
#' pairwise_score(
#'   list("A" = gr_a1, "B" = gr_a2),
#'   list("X" = gr_b1, "Y" = gr_b2),
#'   score_func
#' )
pairwise_score <- function(grlist_1, grlist_2, score_func) {
  if (is.null(names(grlist_1)) || is.null(names(grlist_2))) {
    stop("GRanges lists must be named")
  }
  combinations <- expand.grid(names(grlist_1), names(grlist_2))
  names(combinations) <- c("gr1", "gr2")

  result <- mapply(
    score_func,
    grlist_1[combinations$gr1],
    grlist_2[combinations$gr2],
    USE.NAMES = FALSE,
    SIMPLIFY = TRUE
  )
  cbind(combinations, score = result)
}
