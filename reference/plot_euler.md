# Plot a Euler diagram with intersections between GRanges objects

Further parameters can be passed to the plot function (see plot.euler)

## Usage

``` r
plot_euler(
  grlist,
  names = NULL,
  ignore.strand = TRUE,
  fills = NULL,
  shape = "circle"
)
```

## Arguments

- grlist:

  List of GRanges objects

- names:

  Names to give the sets

- ignore.strand:

  Whether to ignore strand for intersections

- fills:

  Color fills for each set. Eulerr interpolates the intersections

- shape:

  Which shape to use, either circle or ellipse. Circle is the default,
  but ellipse is more likely to be accurate for a difficult fit.

## Value

A euler.diagram object that can be further processed, or plotted via
print, ggsave

## Examples

``` r
gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 18), strand = "+")
plot_euler(list(gr_1, gr_2, gr_3), names = c("a", "b", "c"))
```
