# bedscout 0.0.0.9008 - 02/04/2026

* gr_universe now returns a source column that is a comma separated list of 
the provenance label. If no label is assigned, generic labels will be added.
* impute_feature, annotate_nearest_features, annotate_overlapping_features,
annotate_nearby_features, all return a "annotation" column instead of the old
"nearby_features".


# bedscout 0.0.0.9007 - 01/04/2026

* impute_feature now has a with_ties parameter to refine edge case behavior 
where several loci have the same overlapping score.

# bedscout 0.0.0.9006 - 26/03/2026

* annotate_nearest_features function that reports closest feature to query 
GenomicRanges object and distance.

# bedscout 0.0.0.9005 - 22/07/2024

* vignettes now use default format HTML vignette.
* prettydoc suggestion is now removed.

# bedscout 0.0.0.9004 - 22/07/2024

* peaks_universe function makes use of loci_consensus to obtain a universe of
peaks across conditions and replicates.

# bedscout 0.0.0.9003 - 18/07/2024

* plot_euler now accepts shape allowing for more complex diagrams.
* fit_euler function added to allow inspection of euler objects for diagnostics.

# bedscout 0.0.0.9002 - 18/07/2024

* loci_consensus function added to calculate consensus GRanges over a list.

# bedscout 0.0.0.9001

* NEWS.md file added
* plot_combinations_score now previously sorts the list, to avoid strange
varying shapes in the resulting heatmap.
