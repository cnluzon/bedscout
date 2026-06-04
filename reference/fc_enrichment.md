# Calculate fold change enrichment of a GRanges object's overlap with another GRanges object based on genome size.

The intuitive idea of this metric is how unlikely is to have that
overlap by chance given the size of the genome and the respective sizes
of the ranges. The value can be understood as x times the random chance.

## Usage

``` r
fc_enrichment(gr1, gr2, genome_size, ignore.strand = TRUE)
```

## Arguments

- gr1:

  A GRanges object

- gr2:

  Another GRanges object

- genome_size:

  Genome size

- ignore.strand:

  If FALSE, only matching strand overlaps are counted

## Value

A numeric value for the fc enrichment

## Details

enrichment = (length(intersection) / length(gr2)) / (length(gr1) /
genome_size)

This was originally taken from how ChromHMM calculates enrichment, and
it is a commutative operation: enrichment(gr1, gr2, size) ==
enrichment(gr2, gr1, size)

## Examples

``` r
gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
fc_enrichment(gr_1, gr_2, 30, ignore.strand = TRUE)
#> [1] 1.487603
```
