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
source(paste(scripts, "plotting_parameters.R", sep="")))

# data import
message("data import")
dat_tab <- read.table(paste(data, "CountTable_raw_wheader.txt", sep=""), sep="\t", stringsAsFactors=F, row.names=1, header=T, check.names=F)
design  <- read.table(paste(data, "design_129E.txt", sep=""),                 sep="\t", stringsAsFactors=F,              header=T, check.names=F)

# sort data table accorindg to AGI codes
idx <- order(rownames(dat_tab))
dat_tab <- dat_tab[idx,]

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
Rep       <- factor(design$Rep,       levels=c("1", "2", "3"))
group     <- factor(design$SampleID,  levels=unique(design$SampleID))

model <- model.matrix( ~ 0 + group + Rep)
colnames(model) <- str_replace(colnames(model), "group", "")
colnames(model) <- str_replace(colnames(model), "-", "")

message("Estimating dispersion")
y <- estimateGLMCommonDisp(y, model)
y <- estimateGLMTrendedDisp(y, model)
y <- estimateGLMTagwiseDisp(y, model)

message("GLM fitting")
fit <- glmFit(y, model)

log2cpm <- cpm(y, prior.count=2, log=TRUE)
write.table(log2cpm, file=paste(data, "log2cpm.txt", sep=""), sep="\t", quote=F, col.names=NA, row.names=T)

# LRT
contrasts <- makeContrasts(
                          P_d4     = (   P_350_d4   -     P_MgCl_d4  ),
                          P_d8     = (   P_350_d8   -     P_MgCl_d8  ),
                          P_d12    = (   P_350_d12  -     P_MgCl_d12 ),
                          P_d16    = (   P_350_d16  -     P_MgCl_d16 ),
                          noP_d4   = ( noP_350_d4   -   noP_MgCl_d4  ),
                          noP_d8   = ( noP_350_d8   -   noP_MgCl_d8  ),
                          noP_d12  = ( noP_350_d12  -   noP_MgCl_d12 ),
                          noP_d16  = ( noP_350_d16  -   noP_MgCl_d16 ),
                          levels=model)
contrast_names <- colnames(contrasts)
n <- length(contrast_names)



# DEG selection
message("LRT for DEG identification")
# LRT for each contrasts
LRT.list <- lapply(1:n, function(x) glmLRT(fit, contrast=contrasts[,x]))
names(LRT.list) <- contrast_names

# logFC and PValue tables
logFC_P.list <- lapply(1:n, function(x) {
	table <- LRT.list[[x]]$table[,c(1,4)]
  table$PValue <- p.adjust(table$PValue, method=p.adj.method)
	colnames(table) <- paste(contrast_names[x], colnames(table), sep="_")
	return(table)
	})
logFC_P <- do.call(data.frame, logFC_P.list)
write.table(logFC_P, file=paste(data, "logFC.P.txt", sep=""), sep="\t", quote=F, col.names=NA, row.names=T)

# Significance picking for each tested model
DE.list <- lapply(1:n, function(x) decideTestsDGE(LRT.list[[x]], adjust.method=p.adj.method, p.value=alpha))
names(DE.list) <- contrast_names

# Number of significant differentially abundant OTUs
total    <- sapply(DE.list, function(x) sum(abs(x)))
induced  <- sapply(DE.list, function(x) sum(x ==  1))
reduced  <- sapply(DE.list, function(x) sum(x == -1))
count <- data.frame(total, induced, reduced)
write.table(count, file=paste(stat, "number_of_DEGs.txt", sep=""), quote=F, row.names=T, col.names=NA, sep="\t")

# significance table
DE <- sapply(1:n, function(x) DE.list[[x]][,1])
colnames(DE) <- contrast_names
write.table(DE, file=paste(stat, "GLM_LRT.whole_gene_significance_table.txt", sep=""), sep="\t", quote=T, row.names=T, col.names=NA)







