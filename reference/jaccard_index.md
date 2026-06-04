# Calculate jaccard index between two GenomicRanges objects.

Jaccard index = length(intersection) / length(union)

## Usage

``` r
jaccard_index(gr1, gr2, ignore.strand = TRUE)
```

## Arguments

- gr1:

  GRanges object

- gr2:

  GRanges object

- ignore.strand:

  If FALSE, only matching strand overlaps are counted

## Value

Numeric value representing the jaccard intersection

## Examples

``` r
gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
jaccard_index(gr_1, gr_2, ignore.strand = TRUE)
#> [1] 0.375
```
