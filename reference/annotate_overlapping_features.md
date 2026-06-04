# Annotates a GRanges object with overlapping features in another GRanges object

Distance is measured as distance from end and start depending whether
the target locus is downstream or upstream.

## Usage

``` r
annotate_overlapping_features(
  gr,
  feat_gr,
  name_field,
  minoverlap = 1L,
  ignore.strand = TRUE
)
```

## Arguments

- gr:

  GRanges object to be annotated

- feat_gr:

  GRanges features to be annotated.

- name_field:

  Name of the column that contains the annotation name in feat_gr

- minoverlap:

  Minimum overlapping (used by GenomicRanges findOverlaps)

- ignore.strand:

  Ignore strand when annotating. Default FALSE

## Value

A GRanges object with a nearby_features field in it containing names of
overlapping features in that distance range (separated by commas)

## Examples

``` r

features_gr <- GenomicRanges::GRanges(
  seqnames = c("chr1", "chr1"),
  IRanges::IRanges(c(10,22), c(20,30)),
  strand = c("-", "-"),
  name = c("Feat_A", "Feat_B")
)

gr <- GenomicRanges::GRanges(
  seqnames = c("chr1"),
  IRanges::IRanges(15, 25),
  strand = "+"
)

annotate_overlapping_features(gr, features_gr, "name", ignore.strand = TRUE)
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames    ranges strand |    annotation
#>          <Rle> <IRanges>  <Rle> |   <character>
#>   [1]     chr1     15-25      + | Feat_A,Feat_B
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
