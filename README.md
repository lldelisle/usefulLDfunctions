# usefulLDfunctions
a R package with the functions I use regularly

## Installation
This package depends on GenomicRanges, rtracklayer, stats, grDevices, utils.
If you GenomicRanges and/or rtracklayer are not installed on your R see Installation of dependencies section.

The easiest way to install usefulLDfunctions is using devtools::install_github() from R:
```
if (!"devtools" %in% installed.packages()){
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
```

## Installation of dependencies
As the installation of Bioconductor package depends on the R version you have, I recommend you to load the 2 following functions in your R console (they are part of the package...)
```
rversionAbove <- function(majorT, minorT = 0){
  # Require base
  myMajorn <- as.numeric(R.version$major)
  majorTn <- as.numeric(majorT)
  if (myMajorn > majorTn){
    return(TRUE)
  } else if (myMajorn < majorTn){
    return(TRUE)
  } else {
    minorTn <- as.numeric(minorT)
    myMinorn <- as.numeric(R.version$minor)
    return(myMinorn >= minorTn)
  }
}

safelyLoadAPackageInCRANorBioconductor <-
  function(myPackage, cranRep = "https://stat.ethz.ch/CRAN/"){
  # Require packages base, utils
  # Require function rversionAbove
  # First try to load the package
  if (!suppressWarnings(eval(parse(text = paste0("require(",
                                                 myPackage,
                                                 ", quietly = T)"))))){
    # Download the list of all CRAN packages
    possiblePackages <- utils::available.packages(repos = cranRep)[, "Package"]
    # Test if it is in CRAN
    if (myPackage %in% possiblePackages){
      # Install it specigying the repo
      # to avoid a window to open to choose the repo
      utils::install.packages(myPackage, repos = cranRep)
    } else {
      # If it is not it should be in bioconductor
      if (rversionAbove(3, 5)){
        # With new versions, you need to use BiocManager
        safelyLoadAPackageInCRANorBioconductor("BiocManager")
        # This function is from BiocManager package
        install(myPackage, update = F, ask = F)
      } else {
        # With older versions you need to source biocLite
        # Sometimes you need https and sometimes http
        tryCatch(source("https://bioconductor.org/biocLite.R"),
                 error = function(e){
                     source("http://bioconductor.org/biocLite.R")
                   })
        biocLite(myPackage, suppressUpdates = T,
                 suppressAutoUpdate = T, ask = F)
      }
    }
    # Now that the package is installed you can load it
    eval(parse(text = paste0("require(", myPackage, ", quietly = T)")))
  }
}
```
Then you can load the missing package(s) by:
```
safelyLoadAPackageInCRANorBioconductor("GenomicRanges")
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
```

## Issues
If you have issues, use the Issues in github or send an email to lucille.delisle\@epfl.ch
