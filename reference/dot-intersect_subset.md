# Intersection of a list of GRanges filtering by index. The use of this is to be able to map by powerset

Intersection of a list of GRanges filtering by index. The use of this is
to be able to map by powerset

## Usage

``` r
.intersect_subset(grlist, indices, ignore.strand)
```

## Arguments

- grlist:

  List of GRanges

- indices:

  Array of integers

- ignore.strand:

  Whether to ignore strand in the intersection

## Value

Size of the intersection between grlist\[indices\]
