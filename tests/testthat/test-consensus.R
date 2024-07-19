
test_that("loci_consensus() returns a correct GRanges object", {
  gr_1 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    IRanges::IRanges(c(1, 7, 12), c(4, 9, 14)),
    strand = "-"
  )
  gr_2 <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(5, 8),
    strand = "+"
  )

  overlap <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(1, 12), c(9, 14)),
    strand = "*"
  )

  consensus <- loci_consensus(list(gr_1, gr_2), min_consensus = 1)
  expect_equal(length(consensus), 2)
  expect_equal(consensus, overlap)
})

test_that("loci_consensus() error for min_consensus < 1", {
  gr_1 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    IRanges::IRanges(c(1, 7, 12), c(4, 9, 14)),
    strand = "-"
  )

  expect_error(
    loci_consensus(list(gr_1, gr_2), min_consensus = 0),
    "min_consensus must be a positive number: 0"
  )
})

test_that("loci_consensus() error for min_consensus > length of GRanges list", {
  gr_1 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    IRanges::IRanges(c(1, 7, 12), c(4, 9, 14)),
    strand = "-"
  )

  gr_2 <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(5, 8),
    strand = "+"
  )

  expect_error(
    loci_consensus(list(gr_1, gr_2), min_consensus = 4),
    regexp = "min_consensus must be smaller or equal to size of the list*"
  )
})

test_that("loci_consensus() with single GRanges min_consensus = 1 returns the same GRanges, except strand", {
  gr_1 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    IRanges::IRanges(c(1, 7, 12), c(4, 9, 14)),
    strand = "-"
  )
  overlap <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    IRanges::IRanges(c(1, 7, 12), c(4, 9, 14)),
    strand = "*"
  )
  consensus <- loci_consensus(list(gr_1), min_consensus = 1)
  expect_equal(consensus, overlap)
})

test_that("loci_consensus() filters by coverage", {
  gr_1 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    IRanges::IRanges(c(1, 7, 12), c(4, 9, 14)),
    strand = "-"
  )
  gr_2 <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(5, 8),
    strand = "+"
  )

  overlap <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(7, 8),
    strand = "*"
  )

  consensus <- loci_consensus(list(gr_1, gr_2), min_consensus = 2)
  expect_equal(length(consensus), 1)
  expect_equal(consensus[1, ], overlap)
})

test_that("loci_consensus() works with named list", {
  gr_1 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    IRanges::IRanges(c(1, 7, 12), c(4, 9, 14)),
    strand = "-"
  )
  gr_2 <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(5, 8),
    strand = "+"
  )

  overlap <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(7, 8),
    strand = "*"
  )

  consensus <- loci_consensus(list("A"=gr_1, "B"=gr_2), min_consensus = 2)
  expect_equal(length(consensus), 1)
  expect_equal(consensus[1, ], overlap)
})

test_that("loci_consensus() returns correct result with resize", {
  gr_1 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    IRanges::IRanges(c(1, 7, 12), c(4, 9, 14)),
    strand = "+"
  )
  gr_2 <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(5, 8),
    strand = "+"
  )
  overlap <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(1, 16),
    strand = "*"
  )

  consensus <- loci_consensus(list(gr_1, gr_2), min_consensus = 1, resize = 5, anchor = "start")
  expect_equal(length(consensus), 1)
  expect_equal(consensus, overlap)
})

