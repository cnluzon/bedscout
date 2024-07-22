
#' Total overlapping base pairs between two GRanges objects.
#'
#' @param gr1 GRanges object
#' @param gr2 GRanges object
#' @param ignore.strand Whether to ignore strand when calculating overlap.
#'
#' @return Integer representing the total bp overlap between two sets of loci.
#' @importFrom GenomicRanges intersect width
#' @importClassesFrom GenomicRanges GRanges
#' @importClassesFrom IRanges IRanges
#' @export
#'
#' @examples
#' gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
#' gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
#' total_bp_overlap(gr_1, gr_2, ignore.strand = TRUE)
#' total_bp_overlap(gr_1, gr_2, ignore.strand = FALSE)
total_bp_overlap <- function(gr1, gr2, ignore.strand = TRUE) {
  sum(
    width(
      GenomicRanges::intersect(gr1, gr2, ignore.strand = ignore.strand)
    )
  )
}

#' Calculate jaccard index between two GenomicRanges objects.
#'
#' Jaccard index = length(intersection) / length(union)
#'
#' @param gr1 GRanges object
#' @param gr2 GRanges object
#' @param ignore.strand If FALSE, only matching strand overlaps are counted
#'
#' @return Numeric value representing the jaccard intersection
#' @importFrom GenomicRanges width
#' @export
#'
#' @examples
#' gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
#' gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
#' jaccard_index(gr_1, gr_2, ignore.strand = TRUE)
jaccard_index <- function(gr1, gr2, ignore.strand = TRUE) {
   gr1_len <- sum(GenomicRanges::width(gr1))
   gr2_len <- sum(GenomicRanges::width(gr2))

   bp_ovp <- total_bp_overlap(gr1, gr2, ignore.strand = ignore.strand)
   jaccard <- bp_ovp / (gr1_len + gr2_len - bp_ovp)
   jaccard
}


#' Calculate fold change enrichment of a GRanges object's overlap with another
#' GRanges object based on genome size.
#'
#' The intuitive idea of this metric is how unlikely is to have that overlap by
#' chance given the size of the genome and the respective sizes of the ranges.
#' The value can be understood as x times the random chance.
#'
#' enrichment = (length(intersection) / length(gr2)) / (length(gr1) / genome_size)
#'
#' This was originally taken from how ChromHMM calculates enrichment, and it is
#' a commutative operation: enrichment(gr1, gr2, size) == enrichment(gr2, gr1, size)
#'
#' @param gr1 A GRanges object
#' @param gr2 Another GRanges object
#' @param genome_size Genome size
#' @param ignore.strand If FALSE, only matching strand overlaps are counted
#'
#' @importFrom GenomicRanges width
#' @return A numeric value for the fc enrichment
#' @export
#'
#' @examples
#' gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
#' gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
#' fc_enrichment(gr_1, gr_2, 30, ignore.strand = TRUE)
fc_enrichment <- function(gr1, gr2, genome_size, ignore.strand = TRUE) {
    overlap <- total_bp_overlap(gr1, gr2, ignore.strand = ignore.strand)
    gr1_size <- sum(GenomicRanges::width(gr1))
    gr2_size <- sum(GenomicRanges::width(gr2))
    (overlap / gr2_size) / (gr1_size / genome_size)
}


#' Calculate the number of loci that overlap between two GRanges objects
#' (ignoring length).
#'
#' Note that this is *NOT a commutative operation, since one locus in gr1 could
#' overlap with > 1 locus in gr2, and viceversa. It is what is used generally
#' on venn diagrams, which is why those do not usually add up completely.
#' However, if there is not highly fragmentation difference, results should not
#' be extremely different.
#'
#' @param gr1 A GRanges object
#' @param gr2 A GRanges object
#' @param ignore.strand If FALSE, not matching strands will not be counted.
#' @param minoverlap Minimum overlap in bp to consider this an overlap
#'
#' @importFrom IRanges subsetByOverlaps
#' @return An integer total number of loci in gr1 that overlap with any locus in gr2
#' @export
#'
#' @examples
#' gr1 <- GenomicRanges::GRanges(
#'   seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-"
#' )
#' gr2 <- GenomicRanges::GRanges(
#'   seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+"
#' )
#'
#' loci_overlap(gr1, gr2)
loci_overlap <- function(gr1, gr2, ignore.strand = TRUE, minoverlap = 1L) {
   length(
     IRanges::subsetByOverlaps(gr1, gr2,
        minoverlap = minoverlap,
        ignore.strand = ignore.strand
      )
    )
}
