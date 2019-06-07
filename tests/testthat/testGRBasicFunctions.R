context("GR Basic functions")
library(usefulLDfunctions)
test_that("bedGraphFromGR gives what is expected", {
  gr <- GenomicRanges::GRanges(seqnames = "chr1",
                               ranges = IRanges::IRanges(start = c(1, 11),
                                                         end = c(10, 12)),
                               score = c(20, 30))
  expect_equal(bedGraphFromGR(gr), data.frame(seqnames = c("chr1", "chr1"),
                                              start = c(0, 10),
                                              end = c(10, 12),
                                              score = c(20, 30)))
})
test_that("grFromBedFile gives what is expected", {
  tests_dir <- system.file("tests", package = "usefulLDfunctions")
  test_bed <- file.path(tests_dir, "test6colWithHeader.bed")
  myRanges <- IRanges::IRanges(start = c(127471197, 127472364, 127473531),
                               end = c(127472363, 127473530, 127474697))
  expectedResult <- GenomicRanges::GRanges(seqnames = rep("chr7", 3),
                                           ranges = myRanges,
                                           strand = rep("+", 3))
  S4Vectors::mcols(expectedResult) <- data.frame(name = c("Pos1",
                                                          "Pos2",
                                                          "Pos3"),
                                                 score = rep(0, 3))
  expect_equal(grFromBedFile(test_bed), expectedResult)
})
test_that("grSortedSimplifiedFromNarrowPeak gives what is expected", {
  tests_dir <- system.file("tests", package = "usefulLDfunctions")
  test_bed <- file.path(tests_dir, "test.narrowPeak")
  load(file.path(tests_dir, "simplifiedOrderedTestNarrowPeak.Rdata"))
  expect_equal(grSortedSimplifiedFromNarrowPeak(test_bed), expectedResult)
})
