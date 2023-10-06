test_that("grlist_intersect() returns correct value", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
  gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 18), strand = "+")

  intersection <- grlist_intersect(list(gr_1, gr_2, gr_3), ignore.strand = TRUE)

  expect_equal(length(intersection), 1)
  expect_equal(sum(GenomicRanges::width(intersection)), 3)
})

test_that("grlist_intersect() returns correct value for empty intersection", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
  gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(27, 30), strand = "+")

  intersection <- grlist_intersect(list(gr_1, gr_2, gr_3))
  expect_equal(length(intersection), 0)
})

test_that("grlist_intersect() does not crash for sucessive empty intersections", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
  gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(27, 30), strand = "+")
  gr_4 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(30, 40), strand = "+")

  intersection <- grlist_intersect(list(gr_1, gr_2, gr_3, gr_4))
  expect_equal(length(intersection), 0)
})
