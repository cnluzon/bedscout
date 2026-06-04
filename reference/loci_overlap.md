# Calculate the number of loci that overlap between two GRanges objects (ignoring length).

Note that this is \*NOT a commutative operation, since one locus in gr1
could overlap with \> 1 locus in gr2, and viceversa. It is what is used
generally on venn diagrams, which is why those do not usually add up
completely. However, if there is not highly fragmentation difference,
results should not be extremely different.

## Usage

``` r
loci_overlap(gr1, gr2, ignore.strand = TRUE, minoverlap = 1L)
```

## Arguments

- gr1:

  A GRanges object

- gr2:

  A GRanges object

- ignore.strand:

  If FALSE, not matching strands will not be counted.

- minoverlap:

  Minimum overlap in bp to consider this an overlap

## Value

An integer total number of loci in gr1 that overlap with any locus in
gr2

## Examples

``` r
gr1 <- GenomicRanges::GRanges(
  seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-"
)
gr2 <- GenomicRanges::GRanges(
  seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+"
)

loci_overlap(gr1, gr2)
#> [1] 0
```
