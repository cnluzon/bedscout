% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_combinations_score}
\alias{plot_combinations_score}
\title{Plot triangular pairwise enrichment heatmap for a single GRanges list}
\usage{
plot_combinations_score(grlist, score_func)
}
\arguments{
\item{grlist}{named List of GRanges objects}

\item{score_func}{Any scoring function (must take a pair of GRanges objects)}
}
\value{
A ggplot object
}
\description{
Plot triangular pairwise enrichment heatmap for a single GRanges list
}
\examples{
gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")

gr_a3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")

score_func <- purrr::partial(fc_enrichment, genome_size = 30, ignore.strand = TRUE)

plot_combinations_score(
  list("A" = gr_a1, "B" = gr_a2, "C" = gr_a3),
  score_func
)
}
