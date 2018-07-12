# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ryohei Thomas Nakano
# nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(org.At.tair.db,     quietly=T, warn.conflicts=F)
library(clusterProfiler,    quietly=T, warn.conflicts=F)
library(stringr,            quietly=T, warn.conflicts=F)
library(ggplot2,            quietly=T, warn.conflicts=F)

#diretoreis
data    <- "/your_data_directory/data/"
scripts <- "/your_data_directory/scripts/"
stat    <- "/your_data_directory/statistics/"
fig     <- "/your_data_directory/fig/"
source(paste(scripts, "plotting_parameters.R", sep="")))

# data import
design  <- read.table(paste(data, "design.txt", sep=""), sep="\t", header=T, check.names=F, stringsAsFactors=F)
k       <-  as.numeric(read.table(paste(stat, "AIC_BIC_best_k.txt", sep=""), sep="\t", header=F, check.names=F, stringsAsFactors=F))
cluster <- read.table(paste(stat, "k_means_", k, ".sorted_hclust.txt", sep=""), sep="\t", header=T, check.names=F, stringsAsFactors=F)


# prepare cluster list
clusters <- unique(cluster$Cluster)
cluster_list <- lapply(clusters, function(x) cluster$ID[cluster$Cluster == x] )
names(cluster_list) <- paste("Cluster", clusters, sep="_")


# compare clusters
cg <- compareCluster(geneCluster=cluster_list,
				fun="enrichGO",
				keyType       = "TAIR",
                OrgDb         = org.At.tair.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)


p <- dotplot(cg) +
	theme(
		axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=8),
		axis.text.y=element_text(size=8),
	    panel.background=element_rect(fill="transparent", colour=NA),
	    plot.background=element_rect(fill="transparent", colour=NA))
ggsave(p, file=paste(fig, "enrichGO_dotplot.pdf", sep=""), bg="transparent", width=7.5, height=9)

x_limits <- levels(p$data$Cluster)[c(7,8,9,11)]
GO_limits <- levels(p$data$Description)[14:25]

p <- p + scale_y_discrete(limits=GO_limits) + scale_x_discrete(limits=rev(x_limits))
ggsave(p, file=paste(fig, "enrichGO_dotplot.defense.pdf", sep=""), bg="transparent", width=6.5, height=4)






