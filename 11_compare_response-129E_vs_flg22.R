# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ryohei Thomas Nakano
# nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(stringr,  quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(psych,    quietly=T, warn.conflicts=F)
library(ggplot2,  quietly=T, warn.conflicts=F)

#diretoreis
data    <- "/your_data_directory/data/"
scripts <- "/your_data_directory/scripts/"
stat    <- "/your_data_directory/statistics/"
fig     <- "/your_data_directory/fig/"
source(paste(scripts, "plotting_parameters.R", sep=""))

# data import
logFC.f <-  as.matrix(read.table(paste(data, "logFC.P.flg22.txt", sep=""),                 sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
logFC.r <-  as.matrix(read.table(paste(data, "logFC.P.txt", sep=""),                       sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
DEG     <-            read.table(paste(stat, "DEG.txt", sep=""),                           sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)
k       <- as.numeric(read.table(paste(stat, "AIC_BIC_best_k.txt", sep=""),                sep="\t", header=F,              check.names=F, stringsAsFactors=F))
cluster <-            read.table(paste(stat, "k_means_", k, ".sorted_hclust.txt", sep=""), sep="\t", header=T,              check.names=F, stringsAsFactors=F)

# target
target_cl <- c(1, 15, 13, 6)	# Defense-related genes, see RNAseq_129E analysis results


# subset logFC_P
idx <- str_detect(colnames(logFC.f), "logFC")
logFC.f <- logFC.f[, idx]

idx <- str_detect(colnames(logFC.f), "JA|phr")
logFC.f <- logFC.f[, !idx]

idx <- str_detect(colnames(logFC.r), "logFC")
logFC.r <- logFC.r[, idx]

# merge
share <- intersect(rownames(logFC.r), rownames(logFC.f))
logFC <- cbind(logFC.r[share, ], logFC.f[share, ])


# melt
melt <- melt(logFC)


# define clusters
idx <- melt$Var1 %in% rownames(DEG)
melt_DEG <- melt[idx, ]
melt_DEG$cluster <- "DEGs"

melt_cluster <- melt_DEG
idx <- match(melt_cluster$Var1, cluster$ID)
melt_cluster$cluster <- cluster$Cluster[idx]

idx <- melt_cluster$cluster %in% target_cl
melt_cluster <- melt_cluster[idx,]


# merge
melt <- do.call(rbind, list(melt_DEG, melt_cluster))


# scaling for plotting
melt$colour <- melt$value

max <- quantile(melt$value, .99)
idx <- melt$colour > max
melt$colour[idx] <- max

min <- quantile(melt$value, .01)
idx <- melt$colour < min
melt$colour[idx] <- min


# labeling
melt$Var2 <- str_replace(melt$Var2, "_logFC", "")
melt$Var2 <- factor(melt$Var2, levels=c("P_d4", "P_d8", "P_d12", "P_d16", "noP_d4", "noP_d8", "noP_d12", "noP_d16", "highPflgCol", "lowPflgCol"))
levels(melt$Var2) <- c("129E_4dpi_highP", "129E_8dpi_highP", "129E_12dpi_highP", "129E_16dpi_highP",
					   "129E_4dpi_lowP",  "129E_8dpi_lowP",  "129E_12dpi_lowP",  "129E_16dpi_lowP",
					   "flg22_12dpi_highP", "flg22_12dpi_lowP")




# saturate 
melt_plot <- melt
melt_sat <- melt

idx <- melt_sat$value > 4
melt_sat$value[idx] <- 4

idx <- melt_sat$value < -4
melt_sat$value[idx] <- -4


idx <- str_detect(melt_sat$Var2, "16dpi|flg22")

p <- ggplot(melt_sat[idx,], aes(x=Var2, y=value, colour=colour)) +
	geom_hline(yintercept=0, size=2, colour=c_black) +
	geom_point(position=position_jitterdodge(), alpha=.5) +
	geom_boxplot(data=melt_plot[idx,], fill="white", outlier.shape=NA) +
	facet_grid(. ~ cluster, labeller=label_both, switch="x") +
	scale_colour_gradient2(low=c_cudo_magenta, mid=c_black, high=c_green, midpoint=0, guide=F) +
	scale_y_continuous(breaks=c(-4, -2, 0, 2, 4), labels=c("<-4", "-2", "0", "2", ">4")) +
	coord_cartesian(ylim = c(-4, 4)) +
	theme_RTN +
	theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) +
	labs(x="", y="Folg changes (log2)")
ggsave(p, file=paste(fig, "boxplot_logFC.DEG.16dpi.sat.pdf", sep=""), width=6, height=4.5, bg="transparent")
ggsave(p, file=paste(fig, "boxplot_logFC.DEG.16dpi.sat.png", sep=""), width=6, height=4.5, bg="transparent")




# one way t test
melt$group <- paste(melt$Var2, melt$cluster, sep="_")
group <- unique(melt$group)

sink(paste(stat, "one_way_t_test.logFC.txt", sep=""))
pvalue <- sapply(group, function(x){
	idx <- melt$group == x
	temp <- melt[idx,]
	t <- t.test(temp$value, mu=0)
	print(t)
	return(t$p.value)
})
sink()

adjP <- p.adjust(pvalue, method="bonferroni")
adjP_df <- data.frame(adjP)
adjP_df$sig <- adjP_df$adjP < 0.05
write.table(adjP_df, file=paste(stat, "one_way_t_test.bonferroni.txt", sep=""), col.names=NA, row.names=T, quote=F, sep="\t")



