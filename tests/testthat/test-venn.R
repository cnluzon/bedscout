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

test_that("calculate_venn_intersections yields correct results", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
  gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(22, 30), strand = "+")

  i <- calculate_venn_intersections(list(gr_1, gr_2, gr_3), ignore.strand = TRUE)
  expect_equal(i, c("A"=1, "B"=1, "C"=1, "A&B"=1, "A&C"=0, "B&C"=1, "A&B&C"=0))
})

test_that("calculate_venn_intersections yields correct results using strand", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
  gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(22, 30), strand = "+")

  i <- calculate_venn_intersections(list(gr_1, gr_2, gr_3), ignore.strand = FALSE)
  expect_equal(i, c("A"=1, "B"=1, "C"=1, "A&B"=0, "A&C"=0, "B&C"=1, "A&B&C"=0))
})

test_that("calculate_venn_intersections gives names", {
  gr_1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")
  gr_3 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(22, 30), strand = "+")

  i <- calculate_venn_intersections(list(gr_1, gr_2, gr_3), names = c("a", "b", "c"), ignore.strand = TRUE)
  expect_equal(i, c("a"=1, "b"=1, "c"=1, "a&b"=1, "a&c"=0, "b&c"=1, "a&b&c"=0))
})
