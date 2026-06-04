# Apply a given score function to all pairs in a list of GRanges objects

Apply a given score function to all pairs in a list of GRanges objects

## Usage

``` r
combinations_score(grlist, score_func)
```

## Arguments

- grlist:

  named List of GRanges objects

- score_func:

  Any scoring function (must take a pair of GRanges objects)

## Value

A data frame

## Examples

``` r
gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")
gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")

# Fix the function parameters
score_func <- purrr::partial(
  fc_enrichment,
  genome_size = 30,
  ignore.strand = TRUE
)

combinations_score(
  list("A" = gr_a1, "B" = gr_a2, "X" = gr_b1),
  score_func
)
#>   gr1 gr2 score
#> 1   A   B   0.0
#> 2   A   X   1.8
#> 3   B   X   1.5
```
