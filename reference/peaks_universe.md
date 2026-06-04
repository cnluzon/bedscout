# Peaks universe for a set of conditions and replicates

Peaks universe for a set of conditions and replicates

## Usage

``` r
peaks_universe(
  peak_files,
  min_consensus,
  resize_pre = NULL,
  resize_across = NULL,
  resize_post = NULL,
  anchor = "center"
)
```

## Arguments

- peak_files:

  A list where each element is a list of peak files

- min_consensus:

  Minimum consensus needed to merge peaks from the same condition
  (replicates)

- resize_pre:

  Resize replicate peaks before calling consensus - Default NULL

- resize_across:

  Resize peaks before merging across conditions. Resizing to a small
  size here will make sure that proximal peaks will not get merged, but
  might incur in certain level of overlap if resize_post is larger.

- resize_post:

  Resize final group of peaks to a fixed size. Default NULL

- anchor:

  Where to anchor for resizing. Defaults to center.

## Value

A GRanges objects with the peaks universe
