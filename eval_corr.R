# evaluate correlation coefficient with normalized phenotype
##################################################
eval_corr <- function(x) {   # x is phenotype
	# normalize “pheno” before estimate the correlation coefficients
	n=ncol(x)
	npheno <- NULL
	for(j in 1:n) {
		y <- (x[, j] - mean(x[, j]))/sd(x[, j])
		npheno <- cbind(npheno, y)
	}
	# estimate correlation coefficients between samples 
	corr <- cor(t(npheno))  # matrix of column-wise correlation coefficients
	return(corr) # return correlation coefficients
}
############################################
#(eg) corr <- eval_corr(pheno)


