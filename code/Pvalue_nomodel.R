#################################################################
# Excutable R-script for evaluating p-values of the correlation test
# with independent and dependent variables by using
# rcorr for Pearson and Spearman, and cor.test for Kendall statistics
#################################################################
# include packages
library(Hmisc)   # rcorr() 
library(MASS)
library(clusterGeneration)   # for rcorrmatrix() 
library(lmf)   	# package "lmf" for nearPD() for covariance matrix

# generate input: correlation coefficient, genotypes, and random interaction 
source("eval_corr.R") # estimate correlation coefficient
args <- commandArgs(TRUE)

# number of samples is the control parameter
n <- as.numeric(args[1])  # number of fixed variable, n=20
p <- as.numeric(args[2])  # number of trials to estimate sample variance, p=20
oracle <- as.numeric(args[3]) # select either random(indep) or dependent case

# mean of p-value vs number of attributes (or sample size)
val <- NULL  
#for(m in seq(10, 500, 10)) {   
for(m in seq(100, 1000, 100)) {   
   # perform 'p' trials to evaluate sample mean and variance
   ans <- NULL
   for(i in 1:p) {
      x <- matrix(, n, m) 
      # select random or dependent case...
      if(oracle==0) {
         # generate random 'x' of n variables and m attributes
         for(j in 1:m) { x[,j] <- rnorm(n) }
      } else {
         # generate dependent 'x' of n variables and m attributes
      	 # generate n*n random correlation matrix
      	 ncorr <- rcorrmatrix(n, alphad=1)
      	 # making positive definite & Cholesky decompo.
      	 nchol <- chol(nearPD(ncorr)) 
      	 for(j in 1:m) { x[,j] <- nchol%*%rnorm(n) }
      }

      # normalize each column of pheno "x" 
      # npheno is a matrix class
      npheno <- NULL   # normalized phenotypes
      for(j in 1:m) {
         y <- (x[, j] - mean(x[, j]))/sd(x[, j])
         npheno <- cbind(npheno, y)
      }

      ############ test of correlation with a matrix
      # rcorr() computes a matrix of Pearson's r (or Spearman's rho rank 
      # correlation coefficients) for all possible pairs of columns of a matrix
      # It also returns p-values.

      #re <- rcorr(t(npheno), type="spearman")  
      #re <- rcorr(t(npheno), type="pearson")  
      #p_vec <- re$P[lower.tri(re$P)]   # vector of p-values re$P

      ### test of correlation using cor.test() for Kendall's statistics
      p_vec <- NULL
      for(a in 1:(n-1)) {
         for(b in (a+1):n) {
            #pp <- cor.test(npheno[a,], npheno[b,], method="pearson")
      	    #pp <- cor.test(npheno[a,], npheno[b,], method="spearman")
      	    pp <- cor.test(npheno[a,], npheno[b,], method="kendall")
      	    # print(pp$p.value)
      	    p_vec <- c(p_vec, pp$p.value)
         }
      }
      ############### end of test

      # evaluate mean of n(n-1)/2 p-values
      ans <- c(ans, mean(p_vec))  
   }  # end of p trials

   val <- rbind(val, c(m, mean(ans), sd(ans)))
} # end of each sample size n

if(oracle==0) {
   #write.csv(val, "pvalue_random_spearman20.csv")
   #write.csv(val, "pvalue_random_pearson20.csv")
   write.csv(val, "pvalue_random_kendall20.csv")
} else {
   #write.csv(val, "pvalue_corr_spearman20.csv")
   #write.csv(val, "pvalue_corr_pearson20.csv")
   write.csv(val, "pvalue_corr_kendall20.csv")
}

########################################
#Rscript Pvalue.R  20  20  0 &    # for random
#Rscript Pvalue.R  20  20  1 &    # for correlated
#################################################################
