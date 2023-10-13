test_that("impute_feature() imputes best jaccard overlapping feature", {
  features_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(10,22), c(20,30)),
    strand = c("-", "-"),
    name = c("Feat_A", "Feat_B")
  )

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(15, 25),
    strand = "+"
  )

  result <- impute_feature(gr, features_gr, "name", ignore.strand = TRUE)
  expect_equal(result$feature[1], "Feat_A")
})

test_that("impute_feature() returns a GRanges object", {
  features_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(10,22), c(20,30)),
    strand = c("-", "-"),
    name = c("Feat_A", "Feat_B")
  )

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(15, 25),
    strand = "+"
  )

  result <- impute_feature(gr, features_gr, "name", ignore.strand = TRUE)
  expect_equal(class(result)[[1]], "GRanges")
})

test_that("impute_feature() returns a NA if no overlaps are found", {
  features_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(10,22), c(20,30)),
    strand = c("-", "-"),
    name = c("Feat_A", "Feat_B")
  )

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr2"),
    IRanges::IRanges(15, 25),
    strand = "+"
  )

  result <- impute_feature(gr, features_gr, "name", ignore.strand = TRUE)
  expect_true(is.na(result$feature[1]))
})

test_that("impute_feature() throws no warnings if no overlapping seqnames", {
  features_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(10,22), c(20,30)),
    strand = c("-", "-"),
    name = c("Feat_A", "Feat_B")
  )

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr2"),
    IRanges::IRanges(15, 25),
    strand = "+"
  )

  expect_silent(impute_feature(gr, features_gr, "name", ignore.strand = TRUE))
})

test_that("impute_feature() returns as many ranges as consulted if no overlaps", {
  features_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(100,220), c(200,300)),
    strand = c("-", "-"),
    name = c("Feat_A", "Feat_B")
  )

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    IRanges::IRanges(c(10, 22, 40), c(20, 30, 50)),
    strand = c("-", "-", "+")
  )

  result <- impute_feature(gr, features_gr, "name", ignore.strand = TRUE)
  expect_equal(length(result), 3)
})

test_that("impute_feature() returns correct values with ignore.strand = FALSE", {
  features_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(10,22), c(20,30)),
    strand = c("+", "-"),
    name = c("Feat_A", "Feat_B")
  )

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(15, 25),
    strand = "-"
  )

  result <- impute_feature(gr, features_gr, "name", ignore.strand = FALSE)
  expect_equal(result$feature[1], "Feat_B")
})

test_that("annotate_nearby_features() gives correct features", {
  features_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(10,22), c(20,30)),
    strand = c("-", "-"),
    name = c("Feat_A", "Feat_B")
  )

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(15, 25),
    strand = "+"
  )

  result <- annotate_nearby_features(
    gr,
    features_gr,
    "name",
    ignore.strand = TRUE
  )

  expect_equal(result$nearby_features[[1]], "Feat_A,Feat_B")
})


test_that("annotate_nearby_features() does not give far overlaps", {
  features_gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1"),
    IRanges::IRanges(c(100,220), c(200,300)),
    strand = c("-", "-"),
    name = c("Feat_A", "Feat_B")
  )

  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1"),
    IRanges::IRanges(15, 25),
    strand = "+"
  )

  result <- annotate_nearby_features(
    gr,
    features_gr,
    "name",
    distance_cutoff = 20,
    ignore.strand = TRUE
  )

  expect_true(is.na(result$nearby_features[[1]]))
})
