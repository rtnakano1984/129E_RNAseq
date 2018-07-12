# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ryohei Thomas Nakano
# nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
# library(limma,    quietly=T, warn.conflicts=F)
library(edgeR,    quietly=T, warn.conflicts=F)
library(stringr,  quietly=T, warn.conflicts=F)

#diretoreis
data    <- "/your_data_directory/data/"
scripts <- "/your_data_directory/scripts/"
stat    <- "/your_data_directory/statistics/"
fig     <- "/your_data_directory/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))


# data import
message("data import")
dat_tab <- read.table(paste(data, "GSE87336_JFC_counts.txt", sep=""), sep="\t", stringsAsFactors=F, row.names=1, header=T, check.names=F)
design  <- read.table(paste(data, "design_flg22.txt", sep=""),              sep="\t", stringsAsFactors=F, header=T, check.names=F)

# DGEList object
x <- as.matrix(dat_tab)
y <- DGEList(counts=x)

# extract genes with a least 100 counts over all samples AND expressed more than 10 counts at least in two samples 
message("Filtering and normalization")
cpm_count_coef <- mean(y$samples$lib.size/1000000)
idx <- rowSums(cpm(y))>(100/cpm_count_coef)  & rowSums(cpm(y)>(10/cpm_count_coef))>=2
y <- y[idx, , keep.lib.sizes=F]

# calculate normalization factor (TMM-normalization)
y <- calcNormFactors(y)


# Create model
rep <- as.factor(design$rep)

group <- factor(design$SampleID, levels=unique(design$SampleID))
model <- model.matrix( ~ 0 + group + rep)
colnames(model) <- str_replace(colnames(model), "group", "")
colnames(model) <- str_replace(colnames(model), "-", "")

message("Estimating dispersion")
y <- estimateGLMCommonDisp(y, model)
y <- estimateGLMTrendedDisp(y, model)
y <- estimateGLMTagwiseDisp(y, model)

message("GLM fitting")
fit <- glmFit(y, model)

# LRT
contrasts <- makeContrasts(
                          highPflgCol = (Col_HighP_flg22 - Col_HighP_noTreatment), 
                          lowPflgCol  = (Col_LowP_flg22  - Col_LowP_noTreatment),
                          highPJACol  = (Col_HighP_MeJA  - Col_HighP_noTreatment), 
                          lowPJACol   = (Col_LowP_MeJA   - Col_LowP_noTreatment),
                          highPflgphr = (phr1phl1_HighP_flg22 - phr1phl1_HighP_noTreatment), 
                          lowPflgphr  = (phr1phl1_LowP_flg22  - phr1phl1_LowP_noTreatment),
                          highPJAphr  = (phr1phl1_HighP_MeJA  - phr1phl1_HighP_noTreatment), 
                          lowPJAphrr  = (phr1phl1_LowP_MeJA   - phr1phl1_LowP_noTreatment),
                          levels=model)
contrast_names <- colnames(contrasts)


# DEG selection
message("LRT for DEG identification")
# LRT for each contrasts
LRT.list <- lapply(contrast_names, function(x) glmLRT(fit, contrast=contrasts[,which(contrast_names == x)]))
names(LRT.list) <- contrast_names

# logFC and PValue tables
logFC_P.list <- lapply(contrast_names, function(x) {
	table <- LRT.list[[which(contrast_names == x)]]$table[,c(1,4)]
	colnames(table) <- paste(x, colnames(table), sep="_")
	return(table)
	})
logFC_P <- do.call(data.frame, logFC_P.list)
write.table(logFC_P, file=paste(data, "logFC.P.flg22.txt", sep=""), sep="\t", quote=F, col.names=NA, row.names=T)












