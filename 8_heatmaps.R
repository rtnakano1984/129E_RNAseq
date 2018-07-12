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
data    <- "/your_data_directory/data/"
scripts <- "/your_data_directory/scripts/"
stat    <- "/your_data_directory/statistics/"
fig     <- "/your_data_directory/fig/"
source(paste(scripts, "plotting_parameters.R", sep="")))

# data import
design   <-             read.table(paste(data, "sorted_design.txt", sep=""),                        sep="\t", header=T,              check.names=F, stringsAsFactors=F)
DEG      <-   as.matrix(read.table(paste(stat, "DEG.txt", sep=""),                                  sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
Z        <-   as.matrix(read.table(paste(data, "sorted_Z.txt", sep=""),                             sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
log2cpm  <-   as.matrix(read.table(paste(data, "log2cpm.txt", sep=""),                             sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
logFC_P  <-   as.matrix(read.table(paste(data, "logFC.P.txt", sep=""),                             sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
k        <-  as.numeric(read.table(paste(stat, "AIC_BIC_best_k.txt", sep=""),                       sep="\t", header=F,              check.names=F, stringsAsFactors=F))
cluster  <-             read.table(paste(stat, "k_means_", k, ".sorted_hclust.txt", sep=""),        sep="\t", header=T,              check.names=F, stringsAsFactors=F)

# colour scaling for cluster identifier
sorted_clusters <- unique(cluster$Cluster)
colour <- data.frame(cluster=sorted_clusters, colour=NA, stringsAsFactors=F)
max <- max(Z)
min <- min(Z)
idx <- (1:k)%%2 == 0
colour$colour[idx]  <- max
colour$colour[!idx] <- min
idx <- match(cluster$Cluster, colour$cluster)
cluster_melt <- data.frame(Var1=cluster$ID, Var2="Cluster", value=colour$colour[idx], group="Cluster")

# melt Z
melt <- melt(Z)
melt$Var2 <- factor(melt$Var2, levels=design$Treatment)
levels(melt$Var2) <- design$sample

# merge
melt$group <- "Z"
melt <- rbind(melt, cluster_melt)

# plot whole Z
melt$Var1    <- factor(melt$Var1, levels=rownames(Z))

p <- ggplot(melt, aes(x=Var2, y=Var1, fill=value)) +
	geom_raster() +
	scale_fill_gradient2(low="#0068B7", mid="white", high="#F39800", na.value="gray95") +
	labs(x="", y="", fill="Modified Z Score") +
	facet_grid(. ~ group, space="free", scales="free", drop=T) +
	theme_RTN +
	theme(legend.position="top",
		axis.text=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank(),
		strip.background = element_blank(),
		strip.text.x = element_blank())
ggsave(p, file=paste(fig, "heatmap_Z.kmeans.pdf", sep=""), bg="transparent", width=4, height=8)
ggsave(p, file=paste(fig, "heatmap_Z.kmeans.png", sep=""), bg="transparent", width=4, height=8)


# split logFC_P
idx <- str_detect(colnames(logFC_P), "logFC")
logFC <- logFC_P[, idx]
P     <- logFC_P[, !idx]

# plot logFC
max <- quantile(logFC, .999)
idx <- (logFC > max)
logFC[idx] <- max

min <- quantile(logFC, .001)
idx <- (logFC < min)
logFC[idx] <- min

melt_FC <- melt(logFC)
melt_FC$Var2 <- str_replace(melt_FC$Var2, "_logFC", "")
melt_FC$Var2 <- factor(melt_FC$Var2, levels=c("P_d4", "P_d8", "P_d12", "P_d16", "noP_d4", "noP_d8", "noP_d12", "noP_d16"))
melt_FC$Var1 <- factor(melt_FC$Var1, levels=rownames(Z))

p <- ggplot(melt_FC, aes(x=Var2, y=Var1, fill=value)) +
	geom_raster() +
	scale_fill_gradient2(high=c_dark_green, mid="white", low=c_cudo_magenta, midpoint=0, na.value="gray95") +
	labs(x="", y="", fill="logFC") +
	theme_RTN +
	theme(legend.position="top",
		axis.text=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank()) +
	guides(fill=guide_colourbar(barheight=1, barwidth=5))
ggsave(p, file=paste(fig, "heatmap_logFC.kmeans.pdf", sep=""), bg="transparent", width=2, height=8)
ggsave(p, file=paste(fig, "heatmap_logFC.kmeans.png", sep=""), bg="transparent", width=2, height=8)

# plot logFC
melt_P <- melt(P)
melt_P$Var2 <- str_replace(melt_P$Var2, "_PValue", "")
melt_P$Var2 <- factor(melt_P$Var2, levels=c("P_d4", "P_d8", "P_d12", "P_d16", "noP_d4", "noP_d8", "noP_d12", "noP_d16"))
melt_P$Var1 <- factor(melt_P$Var1, levels=rownames(Z))
melt_P$sig <- as.numeric(melt_P$value < alpha)

p <- ggplot(melt_P, aes(x=Var2, y=Var1, fill=factor(sig))) +
	geom_raster() +
	scale_fill_manual(values=c("0"="white", "1"=c_red)) +
	labs(x="", y="", fill="Significance") +
	theme_RTN +
	theme(legend.position="top",
		legend.key=element_rect(colour=c_black),
		axis.text=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank())
	ggsave(p, file=paste(fig, "heatmap_P.kmeans.pdf", sep=""), bg="transparent", width=2, height=8)
ggsave(p, file=paste(fig, "heatmap_P.kmeans.png", sep=""), bg="transparent", width=2, height=8)


# plot log2cpm
idx <- rownames(log2cpm) %in% rownames(Z)
log2cpm <- log2cpm[idx, ]
melt_CPM <- melt(log2cpm)
melt_CPM$Var2 <- factor(melt_CPM$Var2, levels=design$Treatment)
melt_CPM$Var1 <- factor(melt_CPM$Var1, levels=rownames(Z))
p <- ggplot(melt_CPM, aes(x=Var2, y=Var1, fill=value)) +
	geom_raster() +
	scale_fill_gradient(low="white", high=c_red, na.value="gray95") +
	labs(x="", y="", fill="log2 CPM") +
	theme_RTN +
	theme(legend.position="top",
		axis.text=element_blank(),
		axis.line=element_blank(),
		axis.ticks=element_blank()) +
	guides(fill=guide_colourbar(barheight=1, barwidth=10))
ggsave(p, file=paste(fig, "heatmap_log2cpm.kmeans.pdf", sep=""), bg="transparent", width=4, height=8)
ggsave(p, file=paste(fig, "heatmap_log2cpm.kmeans.png", sep=""), bg="transparent", width=4, height=8)





