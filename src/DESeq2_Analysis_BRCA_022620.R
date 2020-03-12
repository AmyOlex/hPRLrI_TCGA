# Amy Olex
# 2/26/20
# Script to perform DESeq2 DEG analyses on BRCA data.
#
# Part of this code comes from the tutorial at http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#htseq
library("DESeq2")
library("dplyr")

run_contrast <- function(st, outdir, outprefix, control_level) {
  
  #### Time for DESeq Contrasts
  ## st_intermediate_expressed
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = st, directory = data_dir, design= ~ condition)
  
  ## filter out low count rows
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  ## set the control condition to the cohort that is not expressed, so label "0"
  dds$condition <- relevel(dds$condition, ref = control_level)
  
  ## Now lets do the DE Analysis
  dds <- DESeq(dds)
  res <- results(dds)
  res_ordered <- res[order(res$padj),]
  
  
  ## Save the Results
  write.table(res_ordered, file = paste(outdir,outprefix, "_DEGResults.tsv", sep=""), sep="\t", quote = FALSE)
  saveRDS(dds, file = paste(outdir,outprefix, "_DDS.Rds"))
  jpeg(filename = paste(outdir,outprefix, "_MAPlot.jpeg", sep=""))
  plotMA(res)
  dev.off()
  
}



setwd("~/Desktop/CCTR_Git_Repos/CClevenger_TCGA-BRCA/src")

data_dir <- "/Volumes/GoogleDrive/My Drive/Active Collaborations/CClevanger/BRCA_TCGA_DEG_Analyses 022420/BRCA_1222_RNASeq_HTSeq-counts_022620"

sample_metadata <- read.delim(file = "/Volumes/GoogleDrive/My Drive/Active Collaborations/CClevanger/BRCA_TCGA_DEG_Analyses 022420/BRCA_MetadataMaster_031220.tsv", header = TRUE)

sample_metadata_HRPos <- sample_metadata[sample_metadata$hormone.status == "Positive",]
sample_metadata_intermediate <- sample_metadata[sample_metadata$intermediate.expressed == "Yes",]
sample_metadata_HRPos_intermediate <- sample_metadata_HRPos[sample_metadata_HRPos$intermediate.expressed == "Yes",]

file_list <- paste(data_dir, as.character(sample_metadata$count.file.name), sep = "")

## Create tertile contrast columns instead of just the high/low data
# all WT
tertile <- ntile(sample_metadata$wt.expression, 3)
sample_metadata$wt.tertile1v23 <- recode(tertile, "1" = "low", "2" = "high", "3" = "high")
sample_metadata$wt.tertile12v3 <- recode(tertile, "1" = "low", "2" = "low", "3" = "high")
write.table(sample_metadata, file = "/Volumes/GoogleDrive/My Drive/Active Collaborations/CClevanger/BRCA_TCGA_DEG_Analyses 022420/Tertile_Contrasts_WTExpression.txt", sep="\t", quote = FALSE)

# all Intermediate
tertile <- ntile(sample_metadata_intermediate$intermediate.expression, 3)
sample_metadata_intermediate$intermediate.tertile1v23 <- recode(tertile, "1" = "low", "2" = "high", "3" = "high")
sample_metadata_intermediate$intermediate.tertile12v3 <- recode(tertile, "1" = "low", "2" = "low", "3" = "high")
write.table(sample_metadata_intermediate, file = "/Volumes/GoogleDrive/My Drive/Active Collaborations/CClevanger/BRCA_TCGA_DEG_Analyses 022420/Tertile_Contrasts_IntermediateExpression.txt", sep="\t", quote = FALSE)

# all HRPos WT
tertile <- ntile(sample_metadata_HRPos$wt.expression, 3)
sample_metadata_HRPos$wt.HRPos.tertile1v23 <- recode(tertile, "1" = "low", "2" = "high", "3" = "high")
sample_metadata_HRPos$wt.HRPos.tertile12v3 <- recode(tertile, "1" = "low", "2" = "low", "3" = "high")
write.table(sample_metadata_HRPos, file = "/Volumes/GoogleDrive/My Drive/Active Collaborations/CClevanger/BRCA_TCGA_DEG_Analyses 022420/Tertile_Contrasts_WTExpression_HRPos.txt", sep="\t", quote = FALSE)

# all HRPos Intermediate
tertile <- ntile(sample_metadata_HRPos_intermediate$intermediate.expression, 3)
sample_metadata_HRPos_intermediate$intermediate.HRPos.tertile1v23 <- recode(tertile, "1" = "low", "2" = "high", "3" = "high")
sample_metadata_HRPos_intermediate$intermediate.HRPos.tertile12v3 <- recode(tertile, "1" = "low", "2" = "low", "3" = "high")
write.table(sample_metadata_HRPos_intermediate, file = "/Volumes/GoogleDrive/My Drive/Active Collaborations/CClevanger/BRCA_TCGA_DEG_Analyses 022420/Tertile_Contrasts_IntermediateExpression_HRPos.txt", sep="\t", quote = FALSE)

## Create focused sample tables for each contrast
# Compare patients who express the intermediate and those who don't
st_intermediate_expressed <- data.frame(sampleName = sample_metadata$my.uuid, fileName = sample_metadata$count.file.name, condition = sample_metadata$intermediate.expressed)

# Compare patients with high and low expression of the WT form
st_wt.tertile1v23 <- data.frame(sampleName = sample_metadata$my.uuid, fileName = sample_metadata$count.file.name, condition = sample_metadata$wt.tertile1v23)
st_wt.tertile12v3 <- data.frame(sampleName = sample_metadata$my.uuid, fileName = sample_metadata$count.file.name, condition = sample_metadata$wt.tertile12v3)

# Compare patients with high and low expression of the Intermediate form
st_intermediate.tertile1v23 <- data.frame(sampleName = sample_metadata_intermediate$my.uuid, fileName = sample_metadata_intermediate$count.file.name, condition = sample_metadata_intermediate$intermediate.tertile1v23)
st_intermediate.tertile12v3 <- data.frame(sampleName = sample_metadata_intermediate$my.uuid, fileName = sample_metadata_intermediate$count.file.name, condition = sample_metadata_intermediate$intermediate.tertile12v3)

### Now do the same for the Hormone receptor positive data.
# Compare HR Positive patients who express the intermediate and those who don't
st_HRPos.intermediate_expressed <- data.frame(sampleName = sample_metadata_HRPos$my.uuid, fileName = sample_metadata_HRPos$count.file.name, condition = sample_metadata_HRPos$intermediate.expressed)

# Compare HR Positive patients with high and low expression of the WT form
st_HRPos.wt.tertile1v23 <- data.frame(sampleName = sample_metadata_HRPos$my.uuid, fileName = sample_metadata_HRPos$count.file.name, condition = sample_metadata_HRPos$wt.HRPos.tertile1v23)
st_HRPos.wt.tertile12v3 <- data.frame(sampleName = sample_metadata_HRPos$my.uuid, fileName = sample_metadata_HRPos$count.file.name, condition = sample_metadata_HRPos$wt.HRPos.tertile12v3)

# Compare HR Positive patients with high and low expression of the Intermediate form
st_HRPos.intermediate.tertile1v23 <- data.frame(sampleName = sample_metadata_HRPos_intermediate$my.uuid, fileName = sample_metadata_HRPos_intermediate$count.file.name, condition = sample_metadata_HRPos_intermediate$intermediate.HRPos.tertile1v23)
st_HRPos.intermediate.tertile12v3 <- data.frame(sampleName = sample_metadata_HRPos_intermediate$my.uuid, fileName = sample_metadata_HRPos_intermediate$count.file.name, condition = sample_metadata_HRPos_intermediate$intermediate.HRPos.tertile12v3)





#### Ok, now we have 10 sample tables:
##8 st_intermediate_expressed : 1102
##9 st_wt.tertile1v23 : 1102
##10 st_wt.tertile12v3 : 1102

##3 st_intermediate.tertile1v23 : 203
##4 st_intermediate.tertile12v3 : 203

##5 st_HRpos_intermediate_expressed : 646
##6 st_HRPos.wt.tertile1v23 : 646
##7 st_HRPos.wt.tertile12v3 : 646

##1 st_HRPos.intermediate.tertile1v23 : 135
##2 st_HRPos.intermediate.tertile12v3 : 135

## Since these contrasts have been taking a long time to load into R I'm going to do the smallest ones first.


## set up variables for run
outdir <- "/Volumes/GoogleDrive/My Drive/Active Collaborations/CClevanger/BRCA_TCGA_DEG_Analyses 022420/results/"

run_contrast(st = st_HRPos.intermediate.tertile1v23, outdir = outdir, outprefix = "HRPos.intermediate.tertile1v23.ctrlLow", control_level = "low")
run_contrast(st = st_HRPos.intermediate.tertile12v3, outdir = outdir, outprefix = "HRPos.intermediate.tertile12v3.ctrlLow", control_level = "low")

run_contrast(st = st_intermediate.tertile1v23, outdir = outdir, outprefix = "intermediate.tertile1v23.ctrlLow", control_level = "low")
run_contrast(st = st_intermediate.tertile12v3, outdir = outdir, outprefix = "intermediate.tertile12v3.ctrlLow", control_level = "low")

run_contrast(st = st_HRpos_intermediate_expressed, outdir = outdir, outprefix = "HRpos_intermediate_expressed.ctrlNo", control_level = "No")
run_contrast(st = st_HRPos.wt.tertile1v23, outdir = outdir, outprefix = "HRPos.wt.tertile1v23.ctrlLow", control_level = "low")
run_contrast(st = st_HRPos.wt.tertile12v3, outdir = outdir, outprefix = "HRPos.wt.tertile12v3.ctrlLow", control_level = "low")

run_contrast(st = st_intermediate_expressed, outdir = outdir, outprefix = "intermediate_expressed.ctrlNo", control_level = "No")
run_contrast(st = st_wt.tertile1v23, outdir = outdir, outprefix = "wt.tertile1v23.ctrlLow", control_level = "low")
run_contrast(st = st_wt.tertile12v3, outdir = outdir, outprefix = "wt.tertile12v3.ctrlLow", control_level = "low")
