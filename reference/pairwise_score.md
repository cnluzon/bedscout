# Apply a given score function to all pairwise combinations of two named GRanges lists

Apply a given score function to all pairwise combinations of two named
GRanges lists

## Usage

``` r
pairwise_score(grlist_1, grlist_2, score_func)
```

## Arguments

- grlist_1:

  A named List of GRanges objects

- grlist_2:

  A named List of GRanges objects

- score_func:

  Any scoring function (must take a pair of GRanges objects)

## Value

A data frame

## Examples

``` r
gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")
gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")
gr_b2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 25), strand = "+")

# Fix the function parameters
score_func <- purrr::partial(
  fc_enrichment,
  genome_size = 30,
  ignore.strand = FALSE
)
pairwise_score(
  list("A" = gr_a1, "B" = gr_a2),
  list("X" = gr_b1, "Y" = gr_b2),
  score_func
)
#>   gr1 gr2 score
#> 1   A   X   1.8
#> 2   B   X   0.0
#> 3   A   Y   0.0
#> 4   B   Y   3.0
```
