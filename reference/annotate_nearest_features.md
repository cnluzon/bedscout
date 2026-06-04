# Annotate the closest feature

This works similar as annotate nearby features but is

## Usage

``` r
annotate_nearest_features(gr, feat_gr, name_field, ignore.strand = TRUE)
```

## Arguments

- gr:

  Query GRanges object to annotate

- feat_gr:

  Features to annotate it with

- name_field:

  Name of the column in feat_gr to get the name from

- ignore.strand:

  Whether to report loci from any strand or only matching. By default it
  is set to true.

## Value

gr with extra column "nearby_features" and "distance" annotated

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

annotate_nearest_features(gr, features_gr, "name")
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand |    annotation  distance
#>          <Rle> <IRanges>  <Rle> |   <character> <integer>
#>   [1]     chr1     15-25      + | Feat_A,Feat_B         0
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
