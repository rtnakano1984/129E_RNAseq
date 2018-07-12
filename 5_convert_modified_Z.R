# scripts to reproduce the analysis and figures from Garrido-Oter et al., 2018
#
# originally by Ryohei Thomas Nakano
# nakano@mpipz.mpg.de

options(warn=-1)

# clean up
rm(list=ls())

#diretoreis
data    <- "/your_data_directory/data/"

# data import
log2cpm <- read.table(paste(data, "log2cpm.txt", sep=""),  sep="\t", header=T, row.names=1, check.names=F, stringsAsFactors=F)

# transform to matrix
log2cpm <- as.matrix(log2cpm)

# Z-score-like transformation
transZ <- function(x) {
	m   <- median(x)
	mad <- mad(x, constant=1.4826)	# median absolute deviation
	meanAD <- mean(abs(x-mean(x)))		# mean absolute deviation
	if(mad == 0) {
		y <- (x-m)/(1.253314*meanAD)
	} else {
		y <- (x - m)/mad
	}
	return(y)
}

Z <- t(apply(log2cpm, 1, transZ))
write.table(Z, file=paste(data, "mod_Z.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)

# saturate outliers
max <- quantile(Z, .99)
idx <- Z > max
Z[idx] <- max

min <- quantile(Z, .01)
idx <- Z < min
Z[idx] <- min

# export
write.table(Z, file=paste(data, "saturated_mod_Z.txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)

