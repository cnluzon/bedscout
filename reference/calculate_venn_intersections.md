# Calculates all intersections in a list of GRanges object for a venn diagram of an arbitrary number of elements.

Note that this can be computationally demanding (subsets grow
exponentially) and probably pairwise intersections would be of more
interest for larger number of sets, so there is a hard limit on 5 sets.

## Usage

``` r
calculate_venn_intersections(grlist, names = NULL, ignore.strand = TRUE)
```

## Arguments

- grlist:

  List of GRanges objects (maximum 4)

- names:

  identifiers for the GRanges objects.

- ignore.strand:

  If FALSE overlaps will be counted only with matching strand

## Value

A named array for all intersections named as names.

## Details

NOTE: Calculates intersection values as number of loci in the
intersection.

## Examples

``` r
gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 18), strand = "+")
intersection <- calculate_venn_intersections(
  list(gr_1, gr_2, gr_3),
  names = c("a", "b", "c"),
  ignore.strand = TRUE
)
```
