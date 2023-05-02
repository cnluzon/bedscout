
# bedscout

<!-- badges: start -->
<!-- badges: end -->

The goal of `bedscout` is to agglutinate relevant helper functions for
useful comparisons we make between BED files using available libraries
such as `GenomicRanges` and `rtracklayer`. The intention is not to
replace their functionality or reinvent the wheel, but simplify the use
of otherwise more complex and generic functions.

In the future I will also include here some visualization and simple
statistics.

## Installation

You can install the development version of bedscout from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("cnluzon/bedscout")
```

## Example

This is a basic example which calculates the total overlap in basepairs
between two `GenomicRanges` objects:

``` r
library(bedscout)
library(GenomicRanges)

gr_1 <- GRanges(seqnames = c("chr1"), IRanges(10, 20), strand = "-")
gr_2 <- GRanges(seqnames = c("chr1"), IRanges(15, 25), strand = "+")

total_bp_overlap(gr_1, gr_2, ignore.strand = TRUE)
#> [1] 6
```
