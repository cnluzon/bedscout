# Total overlapping base pairs between two GRanges objects.

Total overlapping base pairs between two GRanges objects.

## Usage

``` r
total_bp_overlap(gr1, gr2, ignore.strand = TRUE)
```

## Arguments

- gr1:

  GRanges object

- gr2:

  GRanges object

- ignore.strand:

  Whether to ignore strand when calculating overlap.

## Value

Integer representing the total bp overlap between two sets of loci.

## Examples

``` r
gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
total_bp_overlap(gr_1, gr_2, ignore.strand = TRUE)
#> [1] 6
total_bp_overlap(gr_1, gr_2, ignore.strand = FALSE)
#> [1] 0
```
