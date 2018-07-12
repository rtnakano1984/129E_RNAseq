# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Nina Dombrowski
# nd.microbiota@gmail.com

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(stringr,  quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(pipeR,    quietly=T, warn.conflicts=F)
library(ggplot2,  quietly=T, warn.conflicts=F)

#diretoreis
data    <- "/your_data_directory/data/"
scripts <- "/your_data_directory/scripts/"
stat    <- "/your_data_directory/statistics/"
fig     <- "/your_data_directory/fig/"
source(paste(scripts, "plotting_parameters.R", sep="")))


# data import
message("Data import...")
log2cpm <- read.table(paste(data, "log2cpm.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)
logFC_P <- read.table(paste(data, "logFC.P.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)
DE      <- read.table(paste(stat, "GLM_LRT.whole_gene_significance_table.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)


# extracting DEGs
idx_P  <- rowSums(abs(DE)) != 0
idx_FC <- rowSums(abs(logFC_P[, str_detect(colnames(logFC_P), "logFC")]) > log2(FC_threshold)) != 0

DEG <- intersect(rownames(DE)[idx_P], rownames(logFC_P)[idx_FC])
DEG_table <- cbind(log2cpm[match(DEG, rownames(log2cpm)), ], logFC_P[match(DEG, rownames(logFC_P)), ])
write.table(DEG_table, file=paste(stat, "DEG.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)





