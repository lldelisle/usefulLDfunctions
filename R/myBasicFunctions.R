### USEFUL WHEN LOADING PACKAGE IN SCRIPTS ####

#' Test if the R version used is above the given version
#'
#' @param majorT a string or numerical value of the major version to use as comparison (for example 3)
#' @param minorT a string or numerical value of the minor version to use as comparison (for example 5 or "5.0"). By default it is 0.
#' @return \code{TRUE} if you are using a R version which is at least the one privided and \code{FALSE} if your version is below.
#' @export
rversionAbove <- function(majorT, minorT = 0){
  # Require base
  myMajorn <- as.numeric(R.version$major)
  majorTn <- as.numeric(majorT)
  if (myMajorn > majorTn){
    return(TRUE)
  } else if (myMajorn < majorTn){
    return(FALSE)
  } else {
    minorTn <- as.numeric(minorT)
    myMinorn <- as.numeric(R.version$minor)
    return(myMinorn >= minorTn)
  }
}

#' Load a package and install it if it is not installed using CRAN or Bioconductor.
#'
#' @param myPackage a string with the name of your package
#' @param cranRep the repository you want to use, by default it is the Swiss one.
#' @details It will try to load the package with the function require quietly=T.
#' If the package is not installed, it will check if it is in CRAN.
#' If it is not it will use BiocManager if the version of R is above 3.5.0
#' Or use biocLite if the version is bellow.
#' @import utils
#' @export
#' @examples
#' # With a bioconductor
#' safelyLoadAPackageInCRANorBioconductor("GenomicRanges")
#' # With a CRAN
#' safelyLoadAPackageInCRANorBioconductor("UpSetR")
safelyLoadAPackageInCRANorBioconductor <-
  function(myPackage, cranRep = "https://stat.ethz.ch/CRAN/"){
    # Require packages base, utils
    # Require function rversionAbove
    # First try to load the package
    if (!suppressWarnings(eval(parse(text = paste0("require(",
                                                   myPackage,
                                                   ", quietly = T)"))))){
      # Download the list of all CRAN packages
      possiblePackages <- tryCatch(utils::available.packages(repos = cranRep),
                                   error = function(e){
                                     utils::available.packages()
                                   })[, "Package"]
      # Test if it is in CRAN
      if (myPackage %in% possiblePackages){
        # Install it specigying the repo
        # to avoid a window to open to choose the repo
        utils::install.packages(myPackage, repos = cranRep)
      } else {
        # If it is not it should be in bioconductor
        if (rversionAbove(3, 5)){
          # With new versions, you need to use BiocManager
          safelyLoadAPackageInCRANorBioconductor("BiocManager", cranRep=cranRep)
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

#### READ FUNCTIONS ####

#' Put the content of a tab separated file (with or without header gzip or not) in a dataframe
#' From the first line where the number of fields follow the condition cond
#'
#' @param fn the name of the file (tab delimited file with optionnally headers). Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, getwd().
#' @param cond the condition that the number of columns should follow (for example "==4" or ">=3")
#' @param keepQuote logical whether the quotes should be kept useful for gtf files (default is F)
#' @return The dataframe containing the values of the file \code{fn} (the header is removed).
#' @importFrom utils read.delim
.readFileFromConditionOnNcols <- function(fn, cond, keepQuote=F){
  # Require packages base, utils
  # i will be the first line (excluding commented lines) with data (no header)
  i <- 1
  while (TRUE){
    # header is a data.frame containing the i-th line
    # (excluding commented lines)
    header <- tryCatch(utils::read.delim(gzfile(fn), nrows = 1,
                                         h = F, skip = (i - 1),
                                         comment.char = "#"),
                       error = function(e){
                         NULL
                       })
    if (is.null(header)){
      return(NULL)
    }
    if (!eval(parse(text = paste0("ncol(header)", cond)))){
      # if the number of columns does not fill the condition cond
      # the i-th line (excluding comments) is a header
      i <- i + 1
    } else {
      # if the number of columns fills the condition cond
      # the i-th line (excluding comments) is the first one with data
      break
    }
  }
  # change the quote char if keepQuote
  if (keepQuote){
    quoteChar <- ""
  } else {
    # I use the default of read.delim
    quoteChar <- "\""
  }
  # return the data frame from the i-th line (excluding comments)
  return(utils::read.delim(gzfile(fn), h = F,
                           skip = (i - 1),
                           comment.char = "#",
                           quote = quoteChar))
}

#' Put the content of a bedGraph (with or without header gzip or not) in a dataframe
#'
#' @param fn the name of the bedgraph file (tab delimited file with 4 columns with optionnally headers). Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, getwd().
#' @return The dataframe with 4 columns (chr, start, end, score) containing the values of the file \code{fn} (the header is removed).
#' The bedgraph format is 0-based half open and this remains true for the dataframe.
#' @export
#' @examples
#' tests_dir <- system.file("tests", package="usefulLDfunctions")
#' test_bedgraph <- file.path(tests_dir, "testNoHeader.bedgraph")
#'
#' # Load the bedgraph with no header in a dataframe
#' test_bedgraph_as_df <- readBedGraph(test_bedgraph)
#'
#' test_bedgraph <- file.path(tests_dir, "test1lineHeader.bedgraph")
#'
#' # Load the bedgraph with 1 line header in a dataframe
#' test_bedgraph_as_df <- readBedGraph(test_bedgraph)
#'
#' test_bedgraph <- file.path(tests_dir, "testmultiplelineHeader.bedgraph")
#'
#' # Load the bedgraph with multiple lines header in a dataframe
#' test_bedgraph_as_df <- readBedGraph(test_bedgraph)
readBedGraph <- function(fn){
  # Require packages base
  # Require function .readFileFromConditionOnNcols
  temp_df <- .readFileFromConditionOnNcols(fn, "==4")
  # the name of the columns is adjusted
  colnames(temp_df) <- c("chr", "start", "end", "score")
  return(temp_df)
}

#' Put the content of a bed file (with or without header gzip or not) in a dataframe
#'
#' @param fn the name of the bed file (tab delimited file with at least 3 columns with optionnally headers). Each row of the table appears as one line of the file. If it does not contain an absolute path, the file name is relative to the current working directory, getwd().
#' @return The dataframe with at least 3 columns (chr, start, end, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts) containing the values of the file \code{fn} (the header is removed).
#' The bed format is 0-based half open and this remains true for the dataframe.
#' @export
#' @examples
#' tests_dir <- system.file("tests", package="usefulLDfunctions")
#' test_bed <- file.path(tests_dir, "test3colWithoutHeader.bed")
#'
#' # Load the bed with no header in a dataframe
#' test_bed_as_df <- readBed(test_bed)
#'
#' test_bed <- file.path(tests_dir, "test6colWithHeader.bed")
#'
#' # Load the bed with 6 columns and a header in a dataframe
#' test_bed_as_df <- readBed(test_bed)
#'
#' test_bed <- file.path(tests_dir, "test12colWithHeader.bed.gz")
#'
#' # Load the bed with header with 12 columns gziped in a dataframe
#' test_bed_as_df <- readBed(test_bed)
#'
readBed <- function(fn){
  # Require packages base
  # Require function .readFileFromConditionOnNcols
  temp_df <- .readFileFromConditionOnNcols(fn, ">=3")
  # the name of the columns is adjusted
  colnames(temp_df) <- c("chr", "start", "end", "name", "score",
                         "strand", "thickStart", "thickEnd",
                         "itemRgb", "blockCount", "blockSizes",
                         "blockStarts")[1:ncol(temp_df)]
  return(temp_df)
}


##### CHECK FUNCTIONS #####

#' Check if a variable which should contains a path for a file exists and has a value which corresponds to an existing file
#'
#' @param variableFile a string containing the name of the variable to check
#' @param default the result if the variableFile is not assigned or the file does not exists and isRequired is set to FALSE, default is NA
#' @param isRequired a logical value indicating if the function should stop in case the file does not exists, default is TRUE
#' @return The content of \code{variableFile} if it is the path to an existing file or the default value
#' @details If isRequired is set to \code{TRUE} and the content of \code{variableFile} is not valid, it will stop and execute an error.
#' @export
#' @examples
#' # A file which exists:
#' myFile <- ".Rhistory"
#' myFile <- list.files()[1]
#' myCheckedFile <- checkFile("myFile")
#'
#' # A variable which does not exists:
#' myCheckedFile <- checkFile("myImaginaryVariableWhichShouldNotExists",isRequired=FALSE)
#'
#' # A variable which does not exists but with default:
#' myCheckedFile <- checkFile("myImaginaryVariableWhichShouldNotExists",
#'                            default=".Rhistory", isRequired=FALSE)
#'
#' # A variable which does not exists but is required:
#' \dontrun{
#' myCheckedFile <- checkFile("myImaginaryVariableWhichShouldNotExists")
#'}
checkFile <- function(variableFile, default = NA, isRequired = T){
  # Require base
  if (exists(variableFile)){
    fn <- eval(parse(text = variableFile))
    if (file.exists(fn)){
      return(fn)
    } else{
      warning("the file specified in", variableFile,
              " is ", fn, " and does not exists.\n")
    }
  }
  if (isRequired){
    stop(paste(variableFile, "is not defined or does not exists but required."))
  } else{
    return(default)
  }
}

#' Check if a variable which should contains a path for a directory exists and has a value which corresponds to an existing folder
#'
#' @param variableDirectory a string containing the name of the variable to check
#' @param default the result if the variableDirectory is not assigned or the directory does not exists and isRequired is set to FALSE, default is NA
#' @param isRequired a logical value indicating if the function should stop in case the folder does not exists, default is TRUE
#' @return The content of \code{variableDirectory} if it is the path to an existing directory or the default value
#' @details If isRequired is set to \code{TRUE} and the content of \code{variableDirectory} is not valid, it will stop and execute an error.
#' @export
#' @examples
#' # A directory which exists:
#' myDirectory <- "./"
#' myCheckedDirectory <- checkDirectory("myDirectory")
#'
#' # A variable which does not exists:
#' myCheckedDirectory <- checkDirectory("myImaginaryVariableWhichShouldNotExists",isRequired=FALSE)
#'
#' # A variable which does not exists but with default:
#' myCheckedDirectory <- checkDirectory("myImaginaryVariableWhichShouldNotExists",
#'                                      default="./", isRequired=FALSE)
#'
#' # A variable which does not exists but is required:
#' \dontrun{
#' myCheckedDirectory <- checkDirectory("myImaginaryVariableWhichShouldNotExists")
#'}
checkDirectory <- function(variableDirectory, default = NA, isRequired = T){
  # Require base
  if (exists(variableDirectory)){
    fn <- eval(parse(text = variableDirectory))
    if (dir.exists(fn)){
      return(fn)
    } else{
      warning("the directory specified in", variableDirectory,
              " is ", fn, " and does not exists.\n")
    }
  }
  if (isRequired){
    stop(paste(variableDirectory,
               "is not defined or does not exists but required."))
  } else{
    return(default)
  }
}

#' Check if a variable which should contains numerical values exists and all values are numeric
#'
#' @param variableN a string containing the name of the variable to check
#' @param default the result if the \code{variableN} is not assigned or some values are not numeric and isRequired is set to FALSE, default is NA
#' @param isRequired a logical value indicating if the function should stop in case the values are not numeric, default is TRUE
#' @return The content of \code{variableN} if they are all numeric or the default value
#' @details If isRequired is set to \code{TRUE} and the content of \code{variableN} is not valid, it will stop and execute an error.
#' @export
#' @examples
#' # all values are numeric:
#' myNumbers <- 1:3
#' myCheckedNumbers <- checkNumericalValues("myNumbers")
#'
#' # some values are not numeric:
#' \dontrun{
#' myNumbers <- list(1,2,"T")
#' myCheckedNumbers <- checkNumericalValues("myNumbers")
#'}
#'
#' # A variable which does not exists:
#' myCheckedNumbers <- checkNumericalValues("myImaginaryVariableWhichShouldNotExists",
#'                                          isRequired = FALSE)
#'
#' # A variable which does not exists but with default:
#' myCheckedNumbers <- checkNumericalValues("myImaginaryVariableWhichShouldNotExists",
#'                                           default = c(1,2,12), isRequired = FALSE)
#'
#' # A variable which does not exists but is required:
#' \dontrun{
#' myCheckedNumbers <- checkNumericalValues("myImaginaryVariableWhichShouldNotExists")
#'}
checkNumericalValues <- function(variableN, default = NA, isRequired = T){
  # Require base
  if (exists(variableN)){
    valN <- eval(parse(text = variableN))
    if (all(is.numeric(valN))){
      return(valN)
    } else{
      warning("the numbers specified in", variableN,
              " is/are ", paste(valN, collapse = ","),
              " and does not contains only numeric values.\n")
    }
  }
  if (isRequired){
    stop(paste(variableN, "is not defined or not a number but required."))
  } else{
    return(default)
  }
}

#' Check if a variable which should contains strings values exists and return values overlapping with the possible
#'
#' @param variableS a string containing the name of the variable to check.
#' @param possible a vector containing the accepted string values or NA if there is no constrain. Default is NA.
#' @param default the result if the \code{variableS} is not assigned or no value overlap with \code{possible} (if \code{possible} is not NA) and isRequired is set to FALSE. Default is NA.
#' @param isRequired a logical value indicating if the function should stop in case no value overlap with \code{possible} (if \code{possible} is not NA). Default is TRUE.
#' @return The content of \code{variableS} if they are all included in \code{possible} (if \code{possible} is not NA)
#' Or the overlap between its content and \code{possible} (if \code{possible} is not NA).
#' Or the default value
#' @details If isRequired is set to \code{TRUE} and the content of \code{variableS} is not valid, it will stop and execute an error.
#' @export
#' @examples
#' myStrings <- c("one", "two", "three")
#'
#' # No restriction
#' myCheckedStrings <- checkStrings("myStrings")
#'
#' # Restriction which match
#' myCheckedStrings <- checkStrings("myStrings", possible = c("one", "two", "three", "four"))
#'
#' # Restriction which match partially
#' myCheckedStrings <- checkStrings("myStrings", possible = c("one", "two"))
#'
#' # Restriction which does not match
#' \dontrun{
#' myCheckedStrings <- checkStrings("myStrings", possible = c("One", "Two"))
#'}
#'
#' # A variable which does not exists:
#' myCheckedStrings <- checkStrings("myImaginaryVariableWhichShouldNotExists", isRequired = FALSE)
checkStrings <- function(variableS, possible = NA,
                         default = NA, isRequired = T){
  # Require base
  if (exists(variableS)){
    valS <- eval(parse(text = variableS))
    if (length(possible) == 1 && is.na(possible)){
      return(valS)
    } else {
      if (all(valS %in% possible)){
        return(valS)
      } else {
        newValS <- intersect(valS, possible)
        rejValS <- setdiff(valS, newValS)
        warning(paste(rejValS, collapse = ","), " were proposed as ",
                variableS, " but are not possible.\n")
        if (length(newValS) > 0){
          return(newValS)
        }
      }
    }
  }
  if (isRequired){
    stop(paste(variableS, "is not defined or not a part of",
               paste(possible, collapse = ","), "but required."))
  } else {
    return(default)
  }
}

#' Check if a variable which should contains a logical value exists and its content is logical
#'
#' @param variableL a string containing the name of the variable to check
#' @param default the result if the \code{variableL} is not assigned or the content are not logical and isRequired is set to FALSE, default is NA
#' @param isRequired a logical value indicating if the function should stop in case the values are not numeric, default is TRUE
#' @return The first value contained in \code{variableL} if it is logical or the default value
#' @details If isRequired is set to \code{TRUE} and the content of \code{variableL} is not valid, it will stop and execute an error.
#' @export
#' @examples
#' # the content is logical:
#' myLogicalValue <- TRUE
#' myCheckedLogicalValue <- checkLogicalValue("myLogicalValue")
#'
#' myLogicalValue <- NA
#' myCheckedLogicalValue <- checkLogicalValue("myLogicalValue")
#'
#' # my value is not logical:
#' myLogicalValue <- "T"
#' \dontrun{
#' myCheckedLogicalValue <- checkLogicalValue("myLogicalValue")
#'}
#'
#' # A variable which does not exists:
#' myCheckedLogicalValue <- checkLogicalValue("myImaginaryVariableWhichShouldNotExists",
#'                                            isRequired = FALSE)
#'
#' # A variable which does not exists but with default:
#' myCheckedLogicalValue <- checkLogicalValue("myImaginaryVariableWhichShouldNotExists",
#'                                            default = TRUE, isRequired = FALSE)
#'
#' # A variable which does not exists but is required:
#' \dontrun{
#' myCheckedLogicalValue <- checkLogicalValue("myImaginaryVariableWhichShouldNotExists")
#' }
checkLogicalValue <- function(variableL, default = NA, isRequired = T){
  # Require base
  if (exists(variableL)){
    valL <- eval(parse(text = variableL))[1]
    if (is.logical(valL)){
      return(valL)
    } else {
      warning("the first value of ", variableL, " is ",
              valL, " and is not a logical value.\n")
    }
  }
  if (isRequired){
    stop(paste(variableL, "is not defined or not logical but required."))
  } else {
    return(default)
  }
}

#' Test if a variable could be used as a color
#'
#' @param colorname a string or a numerical value containing the potential color
#' @return A boolean which would say if the colorname is a valid color name or not.
#' @details A valid color name is a numerical value or a string which is in colors() or a character vector with elements of 7 or 9 characters, "#" followed by the red, blue, green and optionally alpha values in hexadecimal (after rescaling to 0 ... 255).
#' @importFrom grDevices colors
#' @export
isValidColor <- function(colorname){
  if (is.numeric(colorname)){
    return(TRUE)
  }
  if (colorname %in% grDevices::colors()){
    return(TRUE)
  }
  if (is.character(colorname)){
    if (nchar(colorname) == 7 || nchar(colorname) == 9){
      if (substr(colorname, 1, 1) == "#"){
        #I should do other checks
        return(TRUE)
      }
    }
  }
  return(FALSE)
}


#### OTHERS #####

#' Extract multiple elements from a list
#'
#' @param bigList a list from where elements will be extracted
#' @param templNames a vector containing indices or names of the element to extract
#' @return A list containing only the elements named in \code{templNames} or whose indices are in \code{templNames} from \code{bigList}
#' @details The name of the output list will match the name of the \code{bigList} even if the indices were used.
#' @export
#' @examples
#' l<-list('a'=c(1,3),'b'=letters[1:10],'c'="fun")
#' subsetByNamesOrIndices(l,1:2)
#' subsetByNamesOrIndices(l,c("a", "c"))
subsetByNamesOrIndices <- function(bigList, templNames){
  templ <- list()
  for (n in templNames){
    templ <- c(templ, list(bigList[[n]]))
  }
  if (is.integer(templNames)){
    names(templ) <- names(bigList)[templNames]
  } else {
    names(templ) <- templNames
  }
  return(templ)
}

#' Find the longest common substring from the end between 2 words
#'
#' @param word1 first string
#' @param word2 second string
#' @return The longest common substring from the end to \code{word1} and \code{word2} or \code{""} if there is no.
#' @export
#' @examples
#' word1 <- "beautiful"
#' word2 <- "useful"
#' commonEnd(word1, word2)
commonEnd <- function (word1, word2) {
  # s1 are all substrings from the end from length of 1 to whole word1
  s1 <- substring(word1, nchar(word1):1, nchar(word1))
  s2 <- substring(word2, nchar(word2):1, nchar(word2))
  w <- which(s1 %in% s2)
  if(length(w)>0){
    return(s2[max(w)])
  } else {
    return(character(1))
  }
}

#' Find the longest common substring from the start between 2 words
#'
#' @param word1 first string
#' @param word2 second string
#' @return The longest common substring from the start to \code{word1} and \code{word2} or \code{""} if there is no.
#' @export
#' @examples
#' word1 <- "useless"
#' word2 <- "useful"
#' commonStart(word1, word2)
commonStart <- function (word1, word2) {
  # s1 are all substrings from the end from length of 1 to whole word1
  s1 <- substring(word1, 1, 1:nchar(word1))
  s2 <- substring(word2, 1, 1:nchar(word2))
  w <- which(s1 %in% s2)
  if(length(w)>0){
    return(s2[max(w)])
  } else {
    return(character(1))
  }
}


#' Simplify a vector of names removing the common end
#'
#' @param vecOfNames a vector with different names to simplify
#' @return a vector with the same length as \code{vecOfNames} where the end of the names were removed if they were identicals.
#' @export
#' @examples
#' vecOfNames <- c("beautiful", "useful", "painful")
#' simplifiedNamesByEnd(vecOfNames)
simplifiedNamesByEnd <- function(vecOfNames){
  if (length(vecOfNames) <= 1){
    return(vecOfNames)
  }
  curCommonEnd <- vecOfNames[1]
  for (i in 2:length(vecOfNames)){
    curCommonEnd <- commonEnd(vecOfNames[i], curCommonEnd)
  }
  return(gsub(paste0(curCommonEnd,"$"),"",vecOfNames))
}

#' Simplify a vector of names removing the common end
#'
#' @param vecOfNames a vector with different names to simplify
#' @return a vector with the same length as \code{vecOfNames} where the end of the names were removed if they were identicals.
#' @export
#' @examples
#' vecOfNames <- c("old", "oldish", "older")
#' simplifiedNamesByStart(vecOfNames)
simplifiedNamesByStart <- function(vecOfNames){
  if (length(vecOfNames) <= 1){
    return(vecOfNames)
  }
  curCommonStart <- vecOfNames[1]
  for (i in 2:length(vecOfNames)){
    curCommonStart <- commonStart(vecOfNames[i], curCommonStart)
  }
  return(gsub(paste0("^", curCommonStart),"",vecOfNames))
}

#' Simplify a vector of names removing the common start and the common end
#'
#' @param vecOfNames a vector with different names to simplify
#' @return a vector with the same length as \code{vecOfNames} where the start and the end of the names were removed if they were identicals.
#' @export
#' @examples
#' vecOfNames <- c("sampleOne.bed", "sampleTwo.bed", "sampleThree.bed")
#' simplifiedNamesByStartEnd(vecOfNames)
simplifiedNamesByStartEnd <- function(vecOfNames){
  return(simplifiedNamesByStart(simplifiedNamesByEnd(vecOfNames)))
}

#' Show the corner of a matrix
#'
#' @param matrix a matrix or a dataframe
#' @param size the size to display (default is 5)
#' @param corner the corner to display ("topleft", "topright", "bottomleft", or "bottomright"), default is "topleft"
#' @param method the method to select the top or bottom, default is 'minmax', can also be 'headtail'.
#' @return a matrix or a dataframe `size` by `size`.
#' @export
#' @examples
#' myHugeMatrix <- matrix(1:10000, nrow = 100)
#' cornerMat(myHugeMatrix, 10)
#' cornerMat(myHugeMatrix, 3, "bottomleft")
#' mySmallMatrix <- matrix(1:9, nrow = 3)
#' cornerMat(mySmallMatrix)
cornerMat <- function(matrix, size = 5, corner = "topleft", method = "minmax"){
  switch(method,
         "minmax" =  switch(corner,
                            "topleft" = matrix[1:min(nrow(matrix), size), 1:min(ncol(matrix), size)],
                            "bottomleft" = matrix[nrow(matrix) - min(nrow(matrix) - 1, size - 1):0, 1:min(ncol(matrix), size)],
                            "topright" = matrix[1:min(nrow(matrix), size), ncol(matrix) - min(ncol(matrix) - 1, size - 1):0],
                            "bottomright" = matrix[nrow(matrix) - min(nrow(matrix) - 1, size - 1):0, ncol(matrix) - min(ncol(matrix) - 1, size - 1):0],
                            stop("the corner value is not valid.")),
         "headtail" = switch(corner,
                             "topleft" = head(matrix, size)[, head(1:ncol(matrix), size)],
                             "bottomleft" = tail(matrix, size)[, head(1:ncol(matrix), size)],
                             "topright" = head(matrix, size)[, tail(1:ncol(matrix), size)],
                             "bottomright" = tail(matrix, size)[, tail(1:ncol(matrix), size)],
                             stop("the corner value is not valid.")),
         stop("Invalid method. Please choose one of 'minmax', 'headtail'")
  )
}
