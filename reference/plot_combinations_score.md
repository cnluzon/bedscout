# Plot triangular pairwise enrichment heatmap for a single GRanges list

Plot triangular pairwise enrichment heatmap for a single GRanges list

## Usage

``` r
plot_combinations_score(grlist, score_func)
```

## Arguments

- grlist:

  named List of GRanges objects

- score_func:

  Any scoring function (must take a pair of GRanges objects)

## Value

A ggplot object

## Examples

``` r
gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")

gr_a3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")

score_func <- purrr::partial(fc_enrichment, genome_size = 30, ignore.strand = TRUE)

plot_combinations_score(
  list("A" = gr_a1, "B" = gr_a2, "C" = gr_a3),
  score_func
)
```
