##################################################
# evaluate partial correlation coefficient with normalized phenotype
##################################################
eval_pcorr <- function(x) {   # x is phenotype
	# normalize “pheno” before estimate the correlation coefficients
	n=ncol(x)
	npheno <- NULL
	for(j in 1:n) {
		y <- (x[, j] - mean(x[, j]))/sd(x[, j])
		npheno <- cbind(npheno, y)
	}
	# estimate partial correlation coefficients between samples 
	# matrix of column-wise correlation coefficients
	# corr <- pcor(t(npheno), method=“pearson”)  
	corr <- pcor(t(npheno))  
	return(corr) # return correlation coefficients
}
############################################
#(eg) corr <- eval_corr(pheno)
#    pcorr <- corr$estimate
#    write.csv(pcorr, "pcorr.csv")
