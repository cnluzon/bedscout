% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/consensus.R
\name{loci_consensus}
\alias{loci_consensus}
\title{Find a consensus set of loci across a set of GRanges objects.}
\usage{
loci_consensus(grlist, min_consensus = 1, resize = NULL, anchor = "center")
}
\arguments{
\item{grlist}{List of GRanges for loci groups to be compared.}

\item{min_consensus}{If a locus appears in at least min_consensus loci sets,
it will be kept. min_consensus must be a number between 1 and length(grlist)}

\item{resize}{Resize the GRanges objects to a fixed size before the check.
If NULL or 0, the loci are not resized.}

\item{anchor}{If resize, where to anchor. By default is set to center.}
}
\value{
A GRanges object with the consensus list.
}
\description{
Ignores strand.
}
\examples{
gr1 <- GenomicRanges::GRanges(
  seqnames = c("chr1"), IRanges::IRanges(c(11, 20), c(15, 30)), strand = "-"
)
gr2 <- GenomicRanges::GRanges(
  seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+"
)

loci_consensus(list(gr1, gr2), min_consensus = 2)
}
