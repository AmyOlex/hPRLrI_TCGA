## Amy Olex
## 11/25/19
## Salmon and tximport pipeline for BRCA chr5 data/
## Much of this code was copied from http://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html
##
## This script imports the chr5 extracted data for all BRCA samples from the Salmon files. 
## It extracts the read counts and TPM values into a single tab-delim text file.
## Gene names and TCGA identifiers are appended to the TPM and raw count data frames.
## Two data files are written out, one with TPM normalized expression and one with raw read counts.
## 
## This script does nothing else except format the raw expression values.
## It is meant to be run on Hershey in the /home/alolex/src folder.

setwd("/home/alolex/FS8600_Bioinformatics/tcga/BRCA/RNASeq/03_Salmon")

library(tximport)
library(readr)
#library(AnnotationHub)
#library(ensembldb)
#library(RNASeqBits)
#library(NMF)
#library(limma)

## Get all the names and the path of the files.
files_h <- file.path("./quantFiles",list.files("./quantFiles"))

## check to ensure they exist
stopifnot(all(file.exists(files_h)))

## Build a tx2gene object that will summarize all transcript expressions to the gene level
#ah <- AnnotationHub()

#annot_h <- query(ah, patter=c("Homo","EnsDb", "87"))
#EnsDb_h <- annot_h[[1]]
#df_human <- transcripts(EnsDb_h, return.type = "DataFrame")
#tx2gene_h <- df_human[,c("tx_id","gene_id")]

## Import Salmon files into a data frame
#txi_h <- tximport(files_h, type = "salmon", tx2gene = tx2gene_h, ignoreTxVersion = TRUE, importer = function(x) read_tsv(x, col_types="cnnnn"))
txi_h <- tximport(files_h, type = "salmon",importer = function(x) read_tsv(x, col_types="cnnnn"), txOut=TRUE)

## Get raw count data
TPM <- txi_h$abundance
counts <- txi_h$counts
lengths <- txi_h$length

## process the sample names
colnames(TPM) <- sub("_gdc_realn_rehead.chr5.R1.quant.sf", "", sub("./quantFiles/", "", files_h))
colnames(counts) <- sub("_gdc_realn_rehead.chr5.R1.quant.sf", "", sub("./quantFiles/", "", files_h))
colnames(lengths) <- sub("_gdc_realn_rehead.chr5.R1.quant.sf", "", sub("./quantFiles/", "", files_h))

GOI_TPM <- t(TPM[c("ENST00000618457.3","ENST00000619676.3"),])
GOI_counts <- t(counts[c("ENST00000618457.3","ENST00000619676.3"),])
GOI_lengths <- t(lengths[c("ENST00000618457.3","ENST00000619676.3"),])


## Now save the files!
write.table(TPM, file="BRCA_1222_chr5_TPM.txt", quote=FALSE, sep="\t", row.names=TRUE)
write.table(counts, file="BRCA_1222_chr5_rawCounts.txt", quote=FALSE, sep="\t", row.names=TRUE)
write.table(lengths, file="BRCA_1222_chr5_txLengths.txt", quote=FALSE, sep="\t", row.names=TRUE)

write.table(GOI_TPM, file="BRCA_1222_chr5_TPM_GOI.txt", quote=FALSE, sep="\t", row.names=TRUE)
write.table(GOI_counts, file="BRCA_1222_chr5_rawCounts_GOI.txt", quote=FALSE, sep="\t", row.names=TRUE)
write.table(GOI_lengths, file="BRCA_1222_chr5_txLengths_GOI.txt", quote=FALSE, sep="\t", row.names=TRUE)

