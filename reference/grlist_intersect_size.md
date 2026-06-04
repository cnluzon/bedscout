# Calculate the size of the intersection of a list of GRanges objects.

Calculate the size of the intersection of a list of GRanges objects.

## Usage

``` r
grlist_intersect_size(grlist, ignore.strand = TRUE)
```

## Arguments

- grlist:

  A list of GRanges

- ignore.strand:

  If TRUE, overlaps are counted without taking strand into account.

## Examples

``` r
gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 18), strand = "+")
intersection <- grlist_intersect_size(list(gr_1, gr_2, gr_3), ignore.strand = TRUE)
```
