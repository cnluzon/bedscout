# Calculate enrichment combinations of different elements within a list

Calculate enrichment combinations of different elements within a list

## Usage

``` r
combinations_enrichment(grlist, genome_size, ignore.strand = TRUE)
```

## Arguments

- grlist:

  Named list of GRanges objects

- genome_size:

  Genome size

- ignore.strand:

  If true, strand is not considered for overlaps.

## Value

A data frame in long format with the enrichment values calculated per
pair.

## Examples

``` r
gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")
gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")

enrichment <- combinations_enrichment(
  list("A" = gr_a1, "B" = gr_a2, "X" = gr_b1),
  genome_size = 30,
  ignore.strand = FALSE
)
```
