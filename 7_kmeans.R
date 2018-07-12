# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ryohei Thomas Nakano
# nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(stringr,   quietly=T, warn.conflicts=F)
library(dplyr,     quietly=T, warn.conflicts=F)
library(pipeR,     quietly=T, warn.conflicts=F)
library(reshape2,  quietly=T, warn.conflicts=F)
library(ggplot2,   quietly=T, warn.conflicts=F)

#diretoreis

#diretoreis
data    <- "/your_data_directory/data/"
scripts <- "/your_data_directory/scripts/"
stat    <- "/your_data_directory/statistics/"
fig     <- "/your_data_directory/fig/"
source(paste(scripts, "plotting_parameters.R", sep="")))

# functions
repZ <- function(fit, x) {
	ids <- names(fit$cluster[fit$cluster==x])
	
	idx <- rownames(Z) %in% ids
	Z_temp <- Z[idx, ]

	if(length(ids)==1){ medZ <- Z_temp } else { medZ <- apply(Z_temp, 2, median) }

	return(medZ)
}

# data import
design  <-             read.table(paste(data, "sorted_design.txt", sep=""),   sep="\t", header=T,              check.names=F, stringsAsFactors=F)
DEG     <-   as.matrix(read.table(paste(stat, "DEG.txt", sep=""),             sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
Z       <-   as.matrix(read.table(paste(data, "saturated_mod_Z.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
k       <-  as.numeric(read.table(paste(stat, "AIC_BIC_best_k.txt", sep=""),  sep="\t", header=F,              check.names=F, stringsAsFactors=F))

# subset to DEGs
idx <- rownames(Z) %in% rownames(DEG)
Z <- Z[idx, ]

# common melted table of Z
melt <- melt(Z)
melt$Var2 <- factor(melt$Var2, levels=design$Treatment)
levels(melt$Var2) <- design$sample


# k means Clustering using Mfuzz package
fit <- kmeans(Z, centers=k)

# plot individual clusters
repZ_mat <- sapply(1:k, repZ, fit=fit)

# corerlation between clusters
colnames(repZ_mat) <- 1:k
cor <- cor(repZ_mat)

# hclsut of clusters
d <- as.dist(1-cor)
hclust <- hclust(d, "average")
sorted_clusters <- rownames(cor)[hclust$order]

pdf(file=paste(fig, "dendrogram_hclust_of_k_mean_clusters_", k, ".pdf", sep=""))
plot(hclust)
dev.off()

#export sorted Z
idx <- match(melt$Var1, names(fit$cluster))
melt$cluster <- factor(fit$cluster[idx], levels=sorted_clusters)

idx <- order(melt$cluster)
sorted_ids <- unique(melt$Var1[idx])

sorted_Z   <- Z[sorted_ids, ]
write.table(sorted_Z, file=paste(data, "sorted_Z.txt", sep=""), row.names=T, col.names=NA, sep="\t", quote=F)

# export clsuter info
cluster <- data.frame(ID=names(fit$cluster), Cluster=fit$cluster)
cluster$Cluster <- factor(cluster$Cluster, levels=sorted_clusters)
idx <- order(cluster$Cluster)
cluster <- cluster[idx,]
write.table(cluster, file=paste(stat, "k_means_", k, ".sorted_hclust.txt", sep=""), sep="\t", quote=F, col.names=T, row.names=F)

