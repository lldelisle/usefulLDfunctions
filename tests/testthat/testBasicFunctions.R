context("Basic functions")
library(usefulLDfunctions)

# test_that("simple check", {
#   myV <- 1
#   expect_equal(length(ls()),1)
#   expect_equal(exists("myV", where=0),TRUE)
#   expect_equal(eval(parse(text = "myV")),1)
#   expect_equal(checkNumericalValues("myV"),1)
# })

# test_that("checkedDirectory gives what is expected when defined", {
#   myDirectory <- "./"
#   myCheckedDir <- checkDirectory("myDirectory")
#   expect_equal(myCheckedDir, "./")
# })
test_that("checkedDirectory gives what is expected when not defined with no default", {
  expect_equal(checkDirectory("myImaginaryVariableWhichShouldNotExists",
                              isRequired = FALSE), NA)
})
test_that("checkedDirectory gives what is expected when not defined with default", {
  expect_equal(checkDirectory("myImaginaryVariableWhichShouldNotExists",
                              default = "./", isRequired = FALSE),
               "./")
})
test_that("checkedDirectory gives what is expected when not defined with no default but required", {
  expect_error(checkDirectory("myImaginaryVariableWhichShouldNotExists"),
               "myImaginaryVariableWhichShouldNotExists is not defined or does not exists but required.")
})
# test_that("checkedFile gives what is expected when defined", {
#   # tests_dir <- system.file("tests", package="usefulLDfunctions")
#   # myFile <- list.files(tests_dir,full.names = T)[1]
#   myFile <- list.files("../",full.names = T)[1]
#
#   expect_equal(checkFile("myFile"), myFile)
# })
test_that("checkedFile gives what is expected when not defined with no default", {
  expect_equal(checkFile("myImaginaryVariableWhichShouldNotExists",
                         isRequired = FALSE), NA)
})
test_that("checkedFile gives what is expected when not defined with default", {
  expect_equal(checkFile("myImaginaryVariableWhichShouldNotExists",
                         default = ".Rhistory", isRequired = FALSE),
               ".Rhistory")
})
test_that("checkedFile gives what is expected when not defined with no default but required", {
  expect_error(checkFile("myImaginaryVariableWhichShouldNotExists"),
               "myImaginaryVariableWhichShouldNotExists is not defined or does not exists but required.")
})
# test_that("checkedLogicalValue gives what is expected when defined", {
#   myLogicalValue <- TRUE
#
#   expect_equal(checkLogicalValue("myLogicalValue"), TRUE)
# })
# test_that("checkedLogicalValue gives what is expected when defined", {
#   myLogicalValue <- NA
#
#   expect_equal(checkLogicalValue("myLogicalValue"), NA)
# })
# test_that("checkedLogicalValue gives what is expected when defined but not logical", {
#   myLogicalValue <- "TRUE"
#
#   expect_error(checkLogicalValue("myLogicalValue"),"the first value of myLogicalValue is TRUE and is not a logical value.")
# })
test_that("checkedLogicalValue gives what is expected when not defined with no default", {
  expect_equal(checkLogicalValue("myImaginaryVariableWhichShouldNotExists",
                                 isRequired = FALSE), NA)
})
test_that("checkedLogicalValue gives what is expected when not defined with default", {
  expect_equal(checkLogicalValue("myImaginaryVariableWhichShouldNotExists",
                                 default = TRUE, isRequired = FALSE),
               TRUE)
})
test_that("checkedLogicalValue gives what is expected when not defined with no default but required", {
  expect_error(checkLogicalValue("myImaginaryVariableWhichShouldNotExists"),
               "myImaginaryVariableWhichShouldNotExists is not defined or not logical but required.")
})
tests_dir <- system.file("tests", package = "usefulLDfunctions")
test_that("readBed load the bed with no header in a dataframe", {
  test_bed <- file.path(tests_dir, "test3colWithoutHeader.bed")
  expect_equal(readBed(test_bed), data.frame(chr = rep("chr7", 3),
                                             start = c(127471196, 127472363,
                                                       127473530),
                                             end = c(127472363, 127473530,
                                                     127474697)))
})
test_that("readBed load the bed with 6 columns and a header in a dataframe", {
  test_bed <- file.path(tests_dir, "test6colWithHeader.bed")
  expect_equal(readBed(test_bed), data.frame(chr = rep("chr7", 3),
                                             start = c(127471196, 127472363,
                                                       127473530),
                                             end = c(127472363, 127473530,
                                                     127474697),
                                             name = c("Pos1", "Pos2", "Pos3"),
                                             score = rep(0, 3),
                                             strand = rep("+", 3)))
})
test_that("readBed load the bed with 12 columns gziped in a dataframe", {
  test_bed <- file.path(tests_dir, "test12colWithHeader.bed.gz")
  expect_equal(readBed(test_bed),
               data.frame(chr = rep("chr22", 2),
                          start = c(1, 2) * 1e3,
                          end = c(5, 6) * 1e3,
                          name = c("cloneA", "cloneB"),
                          score = c(960, 900),
                          strand = c("+", "-"),
                          thickStart = c(1, 2) * 1e3,
                          thickEnd = c(5, 6) * 1e3,
                          itemRgb = c(0, 0),
                          blockCount = c(2, 2),
                          blockSizes = c("567,488,", "433,399,"),
                          blockStarts = c("0,3512", "0,3601")))
})
test_that("readBedGraph load the bedgraph with no header in a dataframe", {
  test_bedgraph <- file.path(tests_dir, "testNoHeader.bedgraph")
  expect_equal(readBedGraph(test_bedgraph),
               data.frame(chr = rep("chr19", 9),
                          start = 49300000 + c(20, 23, 26, 29, 32, 35,
                                               38, 41, 44) * 100,
                          end = 49300000 + c(23, 26, 29, 32, 35, 38,
                                             41, 44, 47) * 100,
                          score = seq(-1, 1, 0.25)))
})
test_that("readBedGraph load the bedgraph with one line header in a dataframe", {
  test_bedgraph <- file.path(tests_dir, "test1lineHeader.bedgraph")
  expect_equal(readBedGraph(test_bedgraph),
               data.frame(chr = rep("chr19", 9),
                          start = 49300000 + c(20, 23, 26, 29, 32,
                                               35, 38, 41, 44) * 100,
                          end = 49300000 + c(23, 26, 29, 32, 35, 38,
                                             41, 44, 47) * 100,
                          score = seq(-1, 1, 0.25)))
})
test_that("readBedGraph load the bedgraph with multiple lines header in a dataframe", {
  test_bedgraph <- file.path(tests_dir, "testmultiplelineHeader.bedgraph")
  expect_equal(readBedGraph(test_bedgraph),
               data.frame(chr = rep("chr19", 9),
                          start = 49300000 + c(20, 23, 26, 29, 32,
                                               35, 38, 41, 44) * 100,
                          end = 49300000 + c(23, 26, 29, 32, 35, 38,
                                             41, 44, 47) * 100,
                          score = seq(-1, 1, 0.25)))
})
test_that("subsetByNamesOrIndices gives expected results", {
  l <- list("a" = c(1, 3), "b" = letters[1:10], "c" = "fun")
  expect_equal(subsetByNamesOrIndices(l, 1:2),
               list("a" = c(1, 3), "b" = letters[1:10]))
  expect_equal(subsetByNamesOrIndices(l, c("a", "c")),
               list("a" = c(1, 3), "c" = "fun"))
})

test_that("commonEnd gives expected results", {
  word1 <- "beautiful"
  word2 <- "useful"
  expect_equal(commonEnd(word1, word2),"ful")
  expect_equal(commonEnd(paste0(word1,"t"), word2),"")
})


test_that("simplifiedByEnd gives expected results", {
  vecOfNames <- c("beautiful", "useful", "painful")
  expect_equal(simplifiedNamesByEnd(vecOfNames),c("beauti","use","pain"))
  vecOfNamesDoNotMatch <- c("beautifully", "useful", "painful")
  expect_equal(simplifiedNamesByEnd(vecOfNamesDoNotMatch),vecOfNamesDoNotMatch)
  expect_equal(simplifiedNamesByEnd(vecOfNamesDoNotMatch[1]),vecOfNamesDoNotMatch[1])
})

test_that("cornerMat gives expected results", {
  myHugeMatrix <- matrix(1:10000, nrow = 100)
  expect_equal(cornerMat(myHugeMatrix, 2), matrix(c(1, 2, 101, 102), ncol = 2))
  expected <- matrix(rep(c(98, 99, 100), 3) +
                       rep(c(0, 100, 200), each = 3),
                     ncol = 3)
  rownames(expected) <- c(" [98,]", " [99,]", "[100,]")
  expect_equal(cornerMat(myHugeMatrix, 3, "bottomleft"), expected)
  mySmallMatrix <- matrix(1:9, nrow = 3)
  expect_equal(cornerMat(mySmallMatrix), mySmallMatrix)
})
