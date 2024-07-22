
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


test_that("gr_universe() returns correct result in example case", {
# 	  01   5    10   15   20
# 	   ----|----|----|----|
# 	R1    ====   ===
# A	R2  =====          ====
# 	R3    ==          ====
#
# 	R1         ====
# B	R2         ===
# 	R3         ====
#
# All:    ==   ===
# 2/3:    ===  ====    ===
# Any:  =====  ====    ===
# 	  01   5    10   15   20
# 	   ----|----|----|----|

  gr_1 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(4, 11), c(7, 13))
  )
  gr_2 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(2, 17), c(6, 20))
  )
  gr_3 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(4, 16), c(5, 19))
  )
  gr_4 <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(9, 12)
  )
  gr_5 <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(9, 11)
  )
  gr_6 <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(9, 12)
  )

  grlist <- list(
    A = list(
      gr_1, gr_2, gr_3
    ),
    B = list(
      gr_4, gr_5, gr_6
    )
  )

  universe <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(4, 9), c(5, 11))
  )

  universe_23 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    IRanges::IRanges(c(4, 9, 17), c(6, 12, 19))
  )

  universe_any <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    IRanges::IRanges(c(2, 9, 16), c(7, 13, 20))
  )

  computed <- gr_universe(grlist, min_consensus = 3)
  computed_23 <- gr_universe(grlist, min_consensus = 2)
  computed_any <- gr_universe(grlist, min_consensus = 1)

  expect_equal(universe, computed)
  expect_equal(universe_23, computed_23)
  expect_equal(universe_any, computed_any)
})

