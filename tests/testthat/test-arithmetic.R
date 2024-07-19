test_that("total_bp_overlap() skips different strands when ignore.strand = FALSE", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")

  expect_equal(total_bp_overlap(gr_1, gr_2, ignore.strand = FALSE), 0)
})

test_that("total_bp_overlap() ignores strand when ignore.strand = TRUE", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")

  # Note that IRanges are closed intervals:
  # [10, 20] intersection with [15, 25] is [15, 20] which length is 6
  expect_equal(total_bp_overlap(gr_1, gr_2, ignore.strand = TRUE), 6)
})

test_that("jaccard_index() calculates correct value", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")

  expect_equal(jaccard_index(gr_1, gr_2, ignore.strand = TRUE), 0.375)
})

test_that("jaccard_index() calculates correct value when ignore.strand = FALSE", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")

  expect_equal(jaccard_index(gr_1, gr_2, ignore.strand = FALSE), 0)
})

test_that("fc_enrichment() calculates correct value", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
  gsize <- 30

  expect_equal(fc_enrichment(gr_1, gr_2, gsize, ignore.strand = TRUE), 1.487603306)
})

test_that("fc_enrichment() calculates correct value with ignore.strand = FALSE", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
  gsize <- 30

  expect_equal(fc_enrichment(gr_1, gr_2, gsize, ignore.strand = FALSE), 0)
})

test_that("fc_enrichment() increases with higher gsize ignore.strand = FALSE", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
  gsize <- 30

  enrichment_1 <- fc_enrichment(gr_1, gr_2, gsize, ignore.strand = TRUE)
  enrichment_2 <- fc_enrichment(gr_1, gr_2, gsize + 20, ignore.strand = TRUE)
  expect_true(enrichment_1 < enrichment_2)
})

test_that("fc_enrichment() is commutative", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
  gsize <- 30

  enrichment_1 <- fc_enrichment(gr_1, gr_2, gsize, ignore.strand = TRUE)
  enrichment_2 <- fc_enrichment(gr_2, gr_1, gsize, ignore.strand = TRUE)
  expect_equal(enrichment_1, enrichment_2)
})

test_that("loci_overlap() returns correct values", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")

  overlap <- loci_overlap(gr_1, gr_2, ignore.strand = TRUE)
  expect_equal(overlap, 1)
})

test_that("loci_overlap() returns correct values with assymetrical case", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"), IRanges::IRanges(c(15, 18), c(19, 25)), strand = "+")

  overlap_a <- loci_overlap(gr_1, gr_2, ignore.strand = TRUE)
  overlap_b <- loci_overlap(gr_2, gr_1, ignore.strand = TRUE)
  expect_equal(overlap_a, 1)
  expect_equal(overlap_b, 2)
})

test_that("loci_overlap() returns correct values with minoverlap", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"), IRanges::IRanges(c(15, 19), c(18, 25)), strand = "+")

  # This should skip [19-25] overlap with [10-20] as the length of that overlap is 2.
  overlap <- loci_overlap(gr_2, gr_1, ignore.strand = TRUE, minoverlap = 3)
  expect_equal(overlap, 1)
})

test_that("loci_overlap() returns correct values with minoverlap edge case", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1", "chr1"), IRanges::IRanges(c(10, 24), c(20, 30)), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(18, 25), strand = "+")

  # gr_2 [25, 30] overlaps 5bp with gr_1 in TOTAL, but 3, 2, respectively.
  # This is an edge case example to showcase how minoverlap works on GRanges
  # because documentation seems to suggest that it will look at total. Turns
  # out, it does not, so results expected with minoverlap 3, only counts the
  # first overlap.
  overlap <- loci_overlap(gr_1, gr_2, ignore.strand = TRUE, minoverlap = 3)
  expect_equal(overlap, 1)
})
