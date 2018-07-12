# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ryohei Thomas Nakano
# nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

# packages
library(reshape2, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
source("~/Dropbox/scripts/ggplot-themes_RTN.R")

#diretoreis
data    <- "/your_data_directory/data/"
scripts <- "/your_data_directory/scripts/"
stat    <- "/your_data_directory/statistics/"
fig     <- "/your_data_directory/fig/"
source(paste(scripts, "plotting_parameters.R", sep="")))

# functions
kmeansAIC <- function(fit){
	m <- ncol(fit$centers)
	n <- length(fit$cluster)
	k <- nrow(fit$centers)
	D <- fit$tot.withinss
	return(c(
		AIC=(D + 2*m*k),
		BIC=(D + log(n)*m*k))=)
}

# data import
design  <-           read.table(paste(data, "sorted_design.txt", sep=""),   sep="\t", header=T,              check.names=F, stringsAsFactors=F)
DEG     <- as.matrix(read.table(paste(stat, "DEG.txt", sep=""),             sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))
Z       <- as.matrix(read.table(paste(data, "saturated_mod_Z.txt", sep=""), sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F))

# subset to DEGs
idx <- rownames(Z) %in% rownames(DEG)
Z <- Z[idx, ]

# k means clustering with a range of cluster numbers and calculate AIC and BIC
min <- 1
max <- 50
range <- c(min:max)

AIC <- t(sapply(range, function(x) kmeansAIC(kmeans(Z, centers=x)) ))
AIC <- data.frame(n=range, AIC)

# plot AIC/BIC
melt <- melt(AIC, id.vars="n")
p <- ggplot(melt, aes(x=n, y=value, colour=variable)) +
	geom_line() +
	geom_point() +
	labs(colour="", x="Number of clusters", y="") +
	theme_RTN +
	theme(legend.position="top")
ggsave(p, file=paste(fig, "AIC_BIC.pdf", sep=""), width=4.5, height=3.5)

# out put best number of k
k <- range[which.min(AIC$BIC)]
out <- data.frame("#AIC/BIC Likelihodd test",
	paste("#Tested range of clusters:  from ", min, " to ", max, sep=""),
	"#Best k (with lowest BIC) is:      ",
	k)
write.table(out, file=paste(stat, "AIC_BIC_best_k.txt", sep=""), quote=F, row.names=F, col.names=F, sep="\n")

