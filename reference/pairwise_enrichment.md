# Calculate pairwise enrichment matrix between two lists of GRanges objects

Calculate pairwise enrichment matrix between two lists of GRanges
objects

## Usage

``` r
pairwise_enrichment(grlist_1, grlist_2, genome_size, ignore.strand = TRUE)
```

## Arguments

- grlist_1:

  A named list of GRanges

- grlist_2:

  A named list of GRanges

- genome_size:

  Genome size

- ignore.strand:

  Ignore strand to calculate overlaps.

## Value

A data frame in long format with the enrichment values calculated per
pair.

## Examples

``` r
gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")

gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")
gr_b2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 25), strand = "+")

gsize <- 30

enrichment <- pairwise_enrichment(
  list("A" = gr_a1, "B" = gr_a2),
  list("X" = gr_b1, "Y" = gr_b2),
  genome_size = gsize,
  ignore.strand = FALSE
)
```
