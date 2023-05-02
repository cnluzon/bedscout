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
