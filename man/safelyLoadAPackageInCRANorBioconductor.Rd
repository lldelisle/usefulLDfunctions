% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/myBasicFunctions.R
\name{safelyLoadAPackageInCRANorBioconductor}
\alias{safelyLoadAPackageInCRANorBioconductor}
\title{Load a package and install it if it is not installed using CRAN or Bioconductor.}
\usage{
safelyLoadAPackageInCRANorBioconductor(
  myPackage,
  cranRep = "https://stat.ethz.ch/CRAN/"
)
}
\arguments{
\item{myPackage}{a string with the name of your package}

\item{cranRep}{the repository you want to use, by default it is the Swiss one.}
}
\description{
Load a package and install it if it is not installed using CRAN or Bioconductor.
}
\details{
It will try to load the package with the function require quietly=T.
If the package is not installed, it will check if it is in CRAN.
If it is not it will use BiocManager if the version of R is above 3.5.0
Or use biocLite if the version is bellow.
}
\examples{
# With a bioconductor
safelyLoadAPackageInCRANorBioconductor("GenomicRanges")
# With a CRAN
safelyLoadAPackageInCRANorBioconductor("UpSetR")
}
