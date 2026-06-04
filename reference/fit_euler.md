# Find a fit for a set of intersections in the form of an euler object.

Find a fit for a set of intersections in the form of an euler object.

## Usage

``` r
fit_euler(grlist, names = NULL, ignore.strand = TRUE, shape = "circle")
```

## Arguments

- grlist:

  List of GRanges objects

- names:

  Names to give the sets

- ignore.strand:

  Whether to ignore strand for intersections

- shape:

  Which shape to use, either circle or ellipse. Circle is the default,
  but ellipse is more likely to be accurate for a difficult fit.

## Value

An euler object

## Examples

``` r
gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 18), strand = "+")
fit_euler(list(gr_1, gr_2, gr_3), names = c("a", "b", "c"))
#>       original fitted residuals regionError
#> a            0      0         0           0
#> b            0      0         0           0
#> c            0      0         0           0
#> a&b          0      0         0           0
#> a&c          0      0         0           0
#> b&c          0      0         0           0
#> a&b&c        1      1         0           0
#> 
#> diagError: 0 
#> stress:    0 
```
