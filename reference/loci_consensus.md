# Find a consensus set of loci across a set of GRanges objects.

Ignores strand.

## Usage

``` r
loci_consensus(grlist, min_consensus = 1, resize = NULL, anchor = "center")
```

## Arguments

- grlist:

  List of GRanges for loci groups to be compared.

- min_consensus:

  If a locus appears in at least min_consensus loci sets, it will be
  kept. min_consensus must be a number between 1 and length(grlist)

- resize:

  Resize the GRanges objects to a fixed size before the check. If NULL
  or 0, the loci are not resized.

- anchor:

  If resize, where to anchor. By default is set to center.

## Value

A GRanges object with the consensus list.

## Examples

``` r
gr1 <- GenomicRanges::GRanges(
  seqnames = c("chr1"), IRanges::IRanges(c(11, 20), c(15, 30)), strand = "-"
)
gr2 <- GenomicRanges::GRanges(
  seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+"
)

loci_consensus(list(gr1, gr2), min_consensus = 2)
#> GRanges object with 1 range and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1     24-25      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
