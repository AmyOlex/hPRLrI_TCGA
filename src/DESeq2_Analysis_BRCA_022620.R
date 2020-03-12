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

setwd("/Volumes/GoogleDrive/My Drive/Active Collaborations/CClevanger/BRCA_TCGA_DEG_Analyses 022420/")
data_dir <- "/Volumes/GoogleDrive/My Drive/Active Collaborations/CClevanger/BRCA_TCGA_DEG_Analyses 022420/BRCA_1222_RNASeq_HTSeq-counts_022620"
sample_metadata <- read.delim(file = "/Volumes/GoogleDrive/My Drive/Active Collaborations/CClevanger/BRCA_TCGA_DEG_Analyses 022420/BRCA_MetadataMaster_031220.tsv", header = TRUE)
outdir <- "/Volumes/GoogleDrive/My Drive/Active Collaborations/CClevanger/BRCA_TCGA_DEG_Analyses 022420/results/"

#setwd("/home/alolex/data/clients/CClevenger/DESeq_Analysis/")
#data_dir <- "/home/alolex/data/clients/CClevenger/DESeq_Analysis/data"
#sample_metadata <- read.delim(file = "/home/alolex/data/clients/CClevenger/DESeq_Analysis/BRCA_MetadataMaster_031220.tsv", header = TRUE)
#outdir <- "/home/alolex/data/clients/CClevenger/DESeq_Analysis/"


sample_metadata_HRPos <- sample_metadata[sample_metadata$hormone.status == "Positive",]
sample_metadata_intermediate <- sample_metadata[sample_metadata$intermediate.expressed == "Yes",]
sample_metadata_HRPos_intermediate <- sample_metadata_HRPos[sample_metadata_HRPos$intermediate.expressed == "Yes",]

file_list <- paste(data_dir, as.character(sample_metadata$count.file.name), sep = "")

## Create tertile contrast columns instead of just the high/low data
# all WT
tertile <- ntile(sample_metadata$wt.expression, 3)
sample_metadata$wt.tertile1v23 <- recode(tertile, "1" = "low", "2" = "high", "3" = "high")
sample_metadata$wt.tertile12v3 <- recode(tertile, "1" = "low", "2" = "low", "3" = "high")
sample_metadata$wt.tertile1v3 <- recode(tertile, "1" = "low", "2" = "med", "3" = "high")

write.table(sample_metadata, file = paste(outdir, "Tertile_Contrasts_WTExpression.txt", sep=""), sep="\t", quote = FALSE)

# all Intermediate
tertile <- ntile(sample_metadata_intermediate$intermediate.expression, 3)
sample_metadata_intermediate$intermediate.tertile1v23 <- recode(tertile, "1" = "low", "2" = "high", "3" = "high")
sample_metadata_intermediate$intermediate.tertile12v3 <- recode(tertile, "1" = "low", "2" = "low", "3" = "high")
sample_metadata_intermediate$intermediate.tertile1v3 <- recode(tertile, "1" = "low", "2" = "med", "3" = "high")
tertile2 <- ntile(sample_metadata_intermediate$ratio.test, 3)
sample_metadata_intermediate$intermediate.ratio.tertile1v3 <- recode(tertile2, "1" = "low", "2" = "med", "3" = "high")
write.table(sample_metadata_intermediate, file = paste(outdir, "Tertile_Contrasts_IntermediateExpression.txt", sep=""), sep="\t", quote = FALSE)

# all HRPos WT
tertile <- ntile(sample_metadata_HRPos$wt.expression, 3)
sample_metadata_HRPos$wt.HRPos.tertile1v23 <- recode(tertile, "1" = "low", "2" = "high", "3" = "high")
sample_metadata_HRPos$wt.HRPos.tertile12v3 <- recode(tertile, "1" = "low", "2" = "low", "3" = "high")
sample_metadata_HRPos$wt.HRPos.tertile1v3 <- recode(tertile, "1" = "low", "2" = "med", "3" = "high")
write.table(sample_metadata_HRPos, file = paste(outdir, "Tertile_Contrasts_WTExpression_HRPos.txt", sep=""), sep="\t", quote = FALSE)

# all HRPos Intermediate
tertile <- ntile(sample_metadata_HRPos_intermediate$intermediate.expression, 3)
sample_metadata_HRPos_intermediate$intermediate.HRPos.tertile1v23 <- recode(tertile, "1" = "low", "2" = "high", "3" = "high")
sample_metadata_HRPos_intermediate$intermediate.HRPos.tertile12v3 <- recode(tertile, "1" = "low", "2" = "low", "3" = "high")
sample_metadata_HRPos_intermediate$intermediate.HRPos.tertile1v3 <- recode(tertile, "1" = "low", "2" = "med", "3" = "high")
write.table(sample_metadata_HRPos_intermediate, file = paste(outdir, "Tertile_Contrasts_IntermediateExpression_HRPos.txt", sep=""), sep="\t", quote = FALSE)

## Create focused sample tables for each contrast
# Compare patients who express the intermediate and those who don't
st_intermediate.expressed <- data.frame(sampleName = sample_metadata$my.uuid, fileName = sample_metadata$count.file.name, condition = sample_metadata$intermediate.expressed)

# Compare patients with high and low expression of the WT form
st_wt.tertile1v23 <- data.frame(sampleName = sample_metadata$my.uuid, fileName = sample_metadata$count.file.name, condition = sample_metadata$wt.tertile1v23)
st_wt.tertile12v3 <- data.frame(sampleName = sample_metadata$my.uuid, fileName = sample_metadata$count.file.name, condition = sample_metadata$wt.tertile12v3)

sample_metadata_sub <- sample_metadata[sample_metadata$wt.tertile1v3 %in% c("low","high"),]
st_wt.tertile1v3 <- data.frame(sampleName = sample_metadata_sub$my.uuid, fileName = sample_metadata_sub$count.file.name, condition = sample_metadata_sub$wt.tertile1v3)

# Compare patients with high and low expression of the Intermediate form
st_intermediate.tertile1v23 <- data.frame(sampleName = sample_metadata_intermediate$my.uuid, fileName = sample_metadata_intermediate$count.file.name, condition = sample_metadata_intermediate$intermediate.tertile1v23)
st_intermediate.tertile12v3 <- data.frame(sampleName = sample_metadata_intermediate$my.uuid, fileName = sample_metadata_intermediate$count.file.name, condition = sample_metadata_intermediate$intermediate.tertile12v3)

sample_metadata_intermediate_sub <- sample_metadata_intermediate[sample_metadata_intermediate$intermediate.tertile1v3 %in% c("low","high"),]
st_intermediate.tertile1v3 <- data.frame(sampleName = sample_metadata_intermediate_sub$my.uuid, fileName = sample_metadata_intermediate_sub$count.file.name, condition = sample_metadata_intermediate_sub$intermediate.tertile1v3)

sample_metadata_intermediate_ratio_sub <- sample_metadata_intermediate[sample_metadata_intermediate$intermediate.ratio.tertile1v3 %in% c("low","high"),]
st_intermediate.ratio.tertile1v3 <- data.frame(sampleName = sample_metadata_intermediate_ratio_sub$my.uuid, fileName = sample_metadata_intermediate_ratio_sub$count.file.name, condition = sample_metadata_intermediate_ratio_sub$intermediate.ratio.tertile1v3)


### Now do the same for the Hormone receptor positive data.
# Compare HR Positive patients who express the intermediate and those who don't
st_HRPos.intermediate.expressed <- data.frame(sampleName = sample_metadata_HRPos$my.uuid, fileName = sample_metadata_HRPos$count.file.name, condition = sample_metadata_HRPos$intermediate.expressed)

# Compare HR Positive patients with high and low expression of the WT form
st_HRPos.wt.tertile1v23 <- data.frame(sampleName = sample_metadata_HRPos$my.uuid, fileName = sample_metadata_HRPos$count.file.name, condition = sample_metadata_HRPos$wt.HRPos.tertile1v23)
st_HRPos.wt.tertile12v3 <- data.frame(sampleName = sample_metadata_HRPos$my.uuid, fileName = sample_metadata_HRPos$count.file.name, condition = sample_metadata_HRPos$wt.HRPos.tertile12v3)

sample_metadata_HRPos_sub <- sample_metadata_HRPos[sample_metadata_HRPos$wt.HRPos.tertile1v3 %in% c("low","high"),]
st_wt.HRPos.tertile1v3 <- data.frame(sampleName = sample_metadata_HRPos_sub$my.uuid, fileName = sample_metadata_HRPos_sub$count.file.name, condition = sample_metadata_HRPos_sub$wt.HRPos.tertile1v3)


# Compare HR Positive patients with high and low expression of the Intermediate form
st_HRPos.intermediate.tertile1v23 <- data.frame(sampleName = sample_metadata_HRPos_intermediate$my.uuid, fileName = sample_metadata_HRPos_intermediate$count.file.name, condition = sample_metadata_HRPos_intermediate$intermediate.HRPos.tertile1v23)
st_HRPos.intermediate.tertile12v3 <- data.frame(sampleName = sample_metadata_HRPos_intermediate$my.uuid, fileName = sample_metadata_HRPos_intermediate$count.file.name, condition = sample_metadata_HRPos_intermediate$intermediate.HRPos.tertile12v3)

sample_metadata_HRPos_intermediate_sub <- sample_metadata_HRPos_intermediate[sample_metadata_HRPos_intermediate$intermediate.HRPos.tertile1v3 %in% c("low","high"),]
st_HRPos.intermediate.tertile1v3 <- data.frame(sampleName = sample_metadata_HRPos_intermediate_sub$my.uuid, fileName = sample_metadata_HRPos_intermediate_sub$count.file.name, condition = sample_metadata_HRPos_intermediate_sub$intermediate.HRPos.tertile1v3)




#### Ok, now we have 15 sample tables:
## Since these contrasts have been taking a long time to load into R I'm going to do the smallest ones first.

#1  st_HRPos.intermediate.tertile1v23 : 135
print("Running #1  st_HRPos.intermediate.tertile1v23 : 135")
run_contrast(st = st_HRPos.intermediate.tertile1v23, outdir = outdir, outprefix = "HRPos.intermediate.tertile1v23.ctrlLow", control_level = "low")

#2  st_HRPos.intermediate.tertile12v3 : 135
print("Running #2  st_HRPos.intermediate.tertile12v3 : 135")
run_contrast(st = st_HRPos.intermediate.tertile12v3, outdir = outdir, outprefix = "HRPos.intermediate.tertile12v3.ctrlLow", control_level = "low")

#3 st_HRPos.intermediate.tertile1v3 : 90
print("Running #3 st_HRPos.intermediate.tertile1v3 : 90")
run_contrast(st = st_HRPos.intermediate.tertile1v3, outdir = outdir, outprefix = "HRPos.intermediate.tertile1v3.ctrlLow", control_level = "low")

#4  st_intermediate.tertile1v23 : 203
print("Running #4  st_intermediate.tertile1v23 : 203")
run_contrast(st = st_intermediate.tertile1v23, outdir = outdir, outprefix = "intermediate.tertile1v23.ctrlLow", control_level = "low")
#5  st_intermediate.tertile12v3 : 203
print("Running #5  st_intermediate.tertile12v3 : 203")
run_contrast(st = st_intermediate.tertile12v3, outdir = outdir, outprefix = "intermediate.tertile12v3.ctrlLow", control_level = "low")
#6 st_intermediate.tertile1v3 : 135
print("Running #6 st_intermediate.tertile1v3 : 135")
run_contrast(st = st_intermediate.tertile1v3, outdir = outdir, outprefix = "intermediate.tertile1v3.ctrlLow", control_level = "low")

#7 st_intermediate.ratio.tertile1v3 : 135
print("Running #7 st_intermediate.ratio.tertile1v3 : 135")
run_contrast(st = st_intermediate.ratio.tertile1v3, outdir = outdir, outprefix = "intermediate.ratio.tertile1v3.ctrlLow", control_level = "low")

#8  st_HRPos.intermediate.expressed : 646
print("Running #8  st_HRPos.intermediate.expressed : 646")
run_contrast(st = st_HRPos.intermediate.expressed, outdir = outdir, outprefix = "HRpos.intermediate.expressed.ctrlNo", control_level = "No")
#9  st_HRPos.wt.tertile1v23 : 646
print("Running #9  st_HRPos.wt.tertile1v23 : 646")
run_contrast(st = st_HRPos.wt.tertile1v23, outdir = outdir, outprefix = "HRPos.wt.tertile1v23.ctrlLow", control_level = "low")
#10  st_HRPos.wt.tertile12v3 : 646
print("Running #10  st_HRPos.wt.tertile12v3 : 646")
run_contrast(st = st_HRPos.wt.tertile12v3, outdir = outdir, outprefix = "HRPos.wt.tertile12v3.ctrlLow", control_level = "low")
#11 st_wt.HRPos.tertile1v3 : 431
print("Running #11 st_wt.HRPos.tertile1v3 : 431")
run_contrast(st = st_wt.HRPos.tertile1v3, outdir = outdir, outprefix = "HRPos.wt.tertile1v3.ctrlLow", control_level = "low")

#12  st_intermediate.expressed : 1102
print("Running #12  st_intermediate.expressed : 1102")
run_contrast(st = st_intermediate.expressed, outdir = outdir, outprefix = "intermediate.expressed.ctrlNo", control_level = "No")
#13  st_wt.tertile1v23 : 1102
print("Running #13  st_wt.tertile1v23 : 1102")
run_contrast(st = st_wt.tertile1v23, outdir = outdir, outprefix = "wt.tertile1v23.ctrlLow", control_level = "low")
#14  st_wt.tertile12v3 : 1102
print("Running #14  st_wt.tertile12v3 : 1102")
run_contrast(st = st_wt.tertile12v3, outdir = outdir, outprefix = "wt.tertile12v3.ctrlLow", control_level = "low")
#15 st_wt.tertile1v3 : 735
print("Running #15 st_wt.tertile1v3 : 735")
run_contrast(st = st_wt.tertile1v3, outdir = outdir, outprefix = "wt.tertile1v3.ctrlLow", control_level = "low")


