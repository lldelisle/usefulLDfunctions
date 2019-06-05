library(GenomicRanges)
library(rtracklayer)

#' Get a GRange with TSS in UCSC format from a path to a gtf file in Ensembl format
#'
#' @param myGTF a file path to a gtf file with ensembl chromosome annotation (no "chr")
#' @return a Genomic Range with all TSS (entries might be duplicated) in UCSC format (with "chr")
#' @export
#' @examples
#' \dontrun{
#' tss <- getTSSinUCSCFormatFromEnsemblGTF("Homo_sapiens.GRCh38.95.gtf.gz")
#' }
#' tss <- getTSSinUCSCFormatFromEnsemblGTF(paste0(system.file("test", package="usefulLDfunctions"),"/Homo_sapiens.GRCh38.95_491firstLines.gtf.gz"))
getTSSinUCSCFormatFromEnsemblGTF <- function(myGTF){
  # require packages GenomicRanges rtracklayer
  gtf <- readGFF(myGTF, filter = list(type = "exon"))
  # The gtf format is like the GRanges format (1-based closed intervals)
  # The chromosome names in UCSC format begins with chr
  gtf$seqid <- paste0("chr",gtf$seqid)
  # The mitochondrial chr is named MT in ensembl and chrM in UCSC
  gtf$seqid[gtf$seqid == "chrMT"] <- "chrM"
  gtfgr <- makeGRangesFromDataFrame(gtf,keep.extra.columns = T)
  tss <- resize(subset(gtfgr,exon_number=="1"), width=1)
  return(tss)
}
