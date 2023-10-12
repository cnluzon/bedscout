test_that("pairwise_enrichment() crashes if lists do not have names", {
  gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")

  gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(10, 20), strand = "-")
  gr_b2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 25), strand = "+")

  gsize <- 100
  expect_error(
    pairwise_enrichment(list(gr_a1, gr_a2), list(gr_b1, gr_b2), genome_size = gsize),
    "GRanges lists must be named"
  )
})

test_that("pairwise_enrichment() yields correct values", {
  gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
  gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")

  gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")
  gr_b2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 25), strand = "+")

  gsize <- 30

  enrichment <- pairwise_enrichment(
    list("A" = gr_a1, "B" = gr_a2),
    list("X" = gr_b1, "Y" = gr_b2),
    genome_size = gsize,
    ignore.strand = FALSE
  )

  expect_equal(enrichment[enrichment$gr1 == "A" & enrichment$gr2 == "X", "enrichment"], 1.8)
  expect_equal(enrichment[enrichment$gr1 == "A" & enrichment$gr2 == "Y", "enrichment"], 0)
  expect_equal(enrichment[enrichment$gr1 == "B" & enrichment$gr2 == "X", "enrichment"], 0)
  expect_equal(enrichment[enrichment$gr1 == "B" & enrichment$gr2 == "Y", "enrichment"], 3)
})

test_that("combination_enrichment() yields correct values", {
  gr_a1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(11, 20), strand = "-")
  gr_a2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(24, 25), strand = "+")

  gr_b1 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(15, 24), strand = "-")
  gr_b2 <- GenomicRanges::GRanges(seqnames = c("chr1"), IRanges::IRanges(16, 25), strand = "+")

  gsize <- 30

  enrichment <- combinations_enrichment(
    list("A" = gr_a1, "B" = gr_a2, "X" = gr_b1, "Y" = gr_b2),
    genome_size = gsize,
    ignore.strand = FALSE
  )

  expect_equal(enrichment[enrichment$gr1 == "A" & enrichment$gr2 == "X", "enrichment"], 1.8)
  expect_equal(enrichment[enrichment$gr1 == "A" & enrichment$gr2 == "Y", "enrichment"], 0)
  expect_equal(enrichment[enrichment$gr1 == "B" & enrichment$gr2 == "X", "enrichment"], 0)
  expect_equal(enrichment[enrichment$gr1 == "B" & enrichment$gr2 == "Y", "enrichment"], 3)
  expect_equal(enrichment[enrichment$gr1 == "A" & enrichment$gr2 == "B", "enrichment"], 0)
  expect_equal(enrichment[enrichment$gr1 == "X" & enrichment$gr2 == "Y", "enrichment"], 0)
})
