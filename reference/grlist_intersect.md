# Find the intersection of a list of GRanges objects.

Note that if only one element is in the list, it returns the same list,
but reducing overlapping ranges within it to 1. This is to keep
consistency on the venn diagrams, since elements that have overlapping
loci will be counted twice.

## Usage

``` r
grlist_intersect(grlist, ignore.strand = TRUE)
```

## Arguments

- grlist:

  A list of GRanges

- ignore.strand:

  If TRUE, overlaps are counted without taking strand into account.

## Value

A GRanges object containing overlapping loci

## Examples

``` r
gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 18), strand = "+")
intersection <- grlist_intersect(list(gr_1, gr_2, gr_3), ignore.strand = TRUE)
```
