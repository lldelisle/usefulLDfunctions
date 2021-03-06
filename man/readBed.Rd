% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myBasicFunctions.R
\name{readBed}
\alias{readBed}
\title{Put the content of a bed file (with or without header gzip or not) in a dataframe}
\usage{
readBed(fn)
}
\arguments{
\item{fn}{the name of the bed file (tab delimited file with at least 3 columns with optionnally headers). Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, getwd().}
}
\value{
The dataframe with at least 3 columns (chr, start, end, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts) containing the values of the file \code{fn} (the header is removed).
The bed format is 0-based half open and this remains true for the dataframe.
}
\description{
Put the content of a bed file (with or without header gzip or not) in a dataframe
}
\examples{
tests_dir <- system.file("tests", package="usefulLDfunctions")
test_bed <- file.path(tests_dir, "test3colWithoutHeader.bed")

# Load the bed with no header in a dataframe
test_bed_as_df <- readBed(test_bed)

test_bed <- file.path(tests_dir, "test6colWithHeader.bed")

# Load the bed with 6 columns and a header in a dataframe
test_bed_as_df <- readBed(test_bed)

test_bed <- file.path(tests_dir, "test12colWithHeader.bed.gz")

# Load the bed with header with 12 columns gziped in a dataframe
test_bed_as_df <- readBed(test_bed)

}
