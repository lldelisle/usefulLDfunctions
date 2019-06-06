#' Convert a GemomicRanges object (1-based close) in a dataframe like bedgraphs (0-based half open)
#'
#' @param gr a GenomicRanges object
#' @return a data frame with 4 columns (seqnames, start, end, score)
#' @export
#' @examples
#' gr <- GenomicRanges::GRanges(seqnames = "chr1",
#'                              ranges = IRanges::IRanges(start = c(1, 11),
#'                                                        end = c(10, 12)),
#'                              score = c(20, 30))
#' bedGraphFromGR(gr)
#' #   seqnames start end score
#' # 1     chr1     0  10    20
#' # 2     chr1    10  12    30
bedGraphFromGR<-function(gr){
  # require package base
  bdg<-data.frame(gr)
  bdg$start<-bdg$start-1
  return(bdg[,c("seqnames","start","end","score")])
}

#' create a GenomicRanges object from a narrowPeak with only one entry per peak (even if multiple summits) sorted by score
#'
#' @param fn a path to the narrowPeak file (it can be gzipped)
#' @return a Genomic ranges sorted by score decreasing with one line per peak, the summit kept is the one with the higher score
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
grSortedSimplifiedFromNarrowPeak<-function(fn){
  # require package GenomicRanges
  # require function .readFileFromConditionOnNcols
  # I read the file from the first line with 10 columns (I remove the header)
  df <- .readFileFromConditionOnNcols(fn, "== 10")
  colnames(df)<-c("chr","start","end","name","score","strand","signalValue","pValue","qValue","relativeSummit")
  # I want to have only one entry per peak by default you may have multiple summits for one peak
  dfA<-stats::aggregate(x=list(score=df$score),by=list(chr=df$chr,start=df$start,end=df$end),FUN=max)
  dfA<-dfA[order(dfA$score,decreasing = T),]
  myWishedRowNames<-paste0(df$chr,":",df$start,"-",df$end,"_",df$score)
  # You may have multiple summits in the same peak which have the same score
  rownames(df)<-make.unique(myWishedRowNames)
  rownames(dfA)<-paste0(dfA$chr,":",dfA$start,"-",dfA$end,"_",dfA$score)
  # Arbitrary the first one will be used
  dfAwithMoreCols<-df[rownames(dfA),]
  myAwithMoreColsGR<-makeGRangesFromDataFrame(dfAwithMoreCols,keep.extra.columns = T)
  return(myAwithMoreColsGR)
}

#' create a GenomicRanges object from a BED
#'
#' @param fn a path to the BED file (it can be gzipped) (0-based half open)
#' @return a Genomic ranges with intervals (1-based closed)
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
grFromBedFileWithName<-function(fn){
  # require package GenomicRanges
  # require function readBed
  temp.df<-readBed(fn)
  colnames(temp.df)<-c("seqname","start","end","name")
  temp.df$start<-temp.df$start+1
  return(makeGRangesFromDataFrame(temp.df,keep.extra.columns = T))
}

