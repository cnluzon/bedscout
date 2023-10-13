#' Find the intersection of a list of GRanges objects.
#'
#' Note that if only one element is in the list, it returns the same list, but
#' reducing overlapping ranges within it to 1. This is to keep consistency on
#' the venn diagrams, since elements that have overlapping loci will be counted
#' twice.
#'
#' @param grlist A list of GRanges
#' @param ignore.strand If TRUE, overlaps are counted without taking strand into
#'   account.
#'
#' @return A GRanges object containing overlapping loci
#' @importFrom purrr reduce
#' @importFrom GenomicRanges intersect
#' @export
#'
#' @examples
#' gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
#' gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
#' gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 18), strand = "+")
#' intersection <- grlist_intersect(list(gr_1, gr_2, gr_3), ignore.strand = TRUE)
grlist_intersect <- function(grlist, ignore.strand = TRUE) {
  if(length(grlist) == 1) {
    GenomicRanges::reduce(grlist[[1]])
  }
  else {
    purrr::reduce(grlist, GenomicRanges::intersect, ignore.strand = ignore.strand)
  }
}

#' Calculate the size of the intersection of a list of GRanges objects.
#' @inheritParams grlist_intersect
#' @export
#'
#' @examples
#' gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
#' gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
#' gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 18), strand = "+")
#' intersection <- grlist_intersect_size(list(gr_1, gr_2, gr_3), ignore.strand = TRUE)
grlist_intersect_size <- function(grlist, ignore.strand = TRUE) {
  length(grlist_intersect(grlist, ignore.strand = ignore.strand))
}

#' Calculates all intersections in a list of GRanges object for a venn diagram
#' of an arbitrary number of elements.
#'
#' Note that this can be computationally demanding (subsets grow exponentially)
#' and probably pairwise intersections would be of more interest for larger
#' number of sets, so there is a hard limit on 5 sets.
#'
#' NOTE: Calculates intersection values as number of loci in the intersection.
#'
#' @param grlist List of GRanges objects (maximum 4)
#' @param names identifiers for the GRanges objects.
#' @param ignore.strand If FALSE overlaps will be counted only with matching strand
#' @return A named array for all intersections named as names.
#' @export
#'
#' @examples
#' gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
#' gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
#' gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 18), strand = "+")
#' intersection <- calculate_venn_intersections(
#'   list(gr_1, gr_2, gr_3),
#'   names = c("a", "b", "c"),
#'   ignore.strand = TRUE
#' )
calculate_venn_intersections <- function(grlist, names = NULL, ignore.strand = TRUE) {
  if (length(grlist) > 5) {
    stop("Intersect loci only supported for up to 5 GRanges objects.")
  }

    if (is.null(names)) {
    names <- LETTERS[1:length(grlist)]
  }

  combinations <- .powerset(seq(1:length(grlist)))
  result_list <- sapply(
     combinations,
     .intersect_subset,
     grlist = grlist,
     ignore.strand = ignore.strand
  )
  result_names <- unlist(lapply(combinations, .make_venn_names, names = names))
  names(result_list) <- result_names

  result_list
}

## Helpers --------------

#' Calculates all possible subsets of a list of elements
#'
#' @param set List of elements
#'
#' @return A list of vectors
.powerset <- function(set) {
  sets <- lapply(set, function(x) utils::combn(set, x, simplify = F))
  unlist(sets, recursive = F)
}

#' Intersection of a list of GRanges filtering by index.
#' The use of this is to be able to map by powerset
#'
#' @param grlist List of GRanges
#' @param indices Array of integers
#' @param ignore.strand Whether to ignore strand in the intersection
#'
#' @return Size of the intersection between grlist\[indices\]
.intersect_subset <- function(grlist, indices, ignore.strand) {
  grlist_sub <- grlist[indices]
  grlist_intersect_size(grlist_sub, ignore.strand)
}

#' Make names for intersections according to index list
#'
#' @param indices Indices to be used
#' @param names Names list
#'
#' @return names\[indices\] joined by an "&" character
.make_venn_names <- function(indices, names) {
  paste(names[indices], collapse = "&")
}

