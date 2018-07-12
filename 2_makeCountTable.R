# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Nina Dombrowski
# nd.microbiota@gmail.com

options(warn=-1)

# cleanup
rm(list=ls())

# packages
library(stringr, quietly=T, warn.conflicts=F)
library(dplyr,   quietly=T, warn.conflicts=F)
library(pipeR,   quietly=T, warn.conflicts=F)

# directory
map  <- "/your_mappping_output_dir/"
data <- "/your_data_directory/data/"

# load design
design <- read.table(paste(data, "design_129E.txt", sep=""), header=T, sep="\t", stringsAsFactors=F)

# get filenames of all .cov files created with bedtools
samples <- paste("Lib", c(LETTERS, paste("A", LETTERS[1:22], sep="")), sep="")
fnames <- paste(samples, "/accepted_hits.cov", sep="")

count_list <- list()

# for each file get a count per gene
count_list <- lapply(c(1:length(fnames)), function(i) {
	print(i)
	
	# read the file
	tt <- read.table(paste(map, fnames[i], sep=""), sep="\t", as.is=T)

	# extract exons
	tt_e <- tt[tt[,3]=="exon",]

	# extract gene ids
	gids <- tt_e$V9
	gids <- str_replace(gids, "Parent=", "")
	gids <- str_replace(gids, "\\.[0-9]{1,}", "")
	gids <- unlist(gids)
	tt_e$gids <- gids

	sum <- tt_e %>>% group_by(gids) %>>% summarise(counts=sum(V10)) %>>% data.frame
	colnames(sum) <- c("Gene", samples[i])

	idx <- order(sum$gids)
	sum <- sum[idx,]

	return(sum)
})

dat_tab <- do.call(cbind, count_list)

idx <- colnames(dat_tab) != "Gene"
dat_tab <- data.frame(Gene=dat_tab[,1], dat_tab[,idx])

idx <- match(colnames(dat_tab, str_replace(design$labels, "_", ""))
colnames(dat_tab) <- design$Treatment[idx]

write.table(dat_tab, file=paste(data, "CountTable_raw_wheader.txt", sep="") , sep="\t", quote=F, row.names=F, col.names=T)
