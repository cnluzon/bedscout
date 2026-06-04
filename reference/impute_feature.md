# Annotate each locus in a GRanges with the best overlap in another feature GRanges. Best overlap is calculated as jaccard index:

length(intersect) / length(union)

## Usage

``` r
impute_feature(
  gr,
  feature_gr,
  name_field,
  minoverlap = 1L,
  ignore.strand = TRUE,
  with_ties = TRUE
)
```

## Arguments

- gr:

  GRanges to annotate

- feature_gr:

  GRanges with features - Needs to have a name field

- name_field:

  Name of the column that contains the feature names in feature_gr

- minoverlap:

  Minimum overlap to consider.

- ignore.strand:

  Whether to ignore strand on annotation.

- with_ties:

  Return ties or return just one result per range? Defaults to True.

## Value

GRanges gr annotated with a "feature" metadata column with the name of
the feature that best overlaps with it, NA if none do.

## Details

If there are ties in jaccard index, ties_policy decides what is
returned.

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

impute_feature(gr, features_gr, "name", ignore.strand = TRUE)
#> GRanges object with 1 range and 2 metadata columns:
#>       seqnames    ranges strand | impute_score  annotation
#>          <Rle> <IRanges>  <Rle> |    <numeric> <character>
#>   [1]     chr1     15-25      + |        0.375      Feat_A
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
