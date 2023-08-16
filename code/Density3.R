#################################################################
# excutable R-script for evaluating sample variance of graph density and 
# transitivity in correlated and uncorrelated values (or phenotypes) by using
# cor.test (instead of rcorr) for Kendall statistics.
# note that we have generated simulated data without any model
#################################################################
# include packages
library(Hmisc)   # rcorr() 
library(igraph)
library(MASS)
library(clusterGeneration)   # for rcorrmatrix() 
library(lmf)   	# package "lmf" for nearPD() for covariance matrix

# generate input: correlation coefficient, genotypes, and random interaction 
source("eval_corr.R") # estimate correlation coefficient
source("cgraph.R")   # graphical analysis of data

args <- commandArgs(TRUE)

# number of samples is the control parameter
m <- as.numeric(args[1])  # number of phenotypes (attributes)
p <- as.numeric(args[2])  # number of trials to estimate sample variance
oracle <- as.numeric(args[3]) # select either random or correlated case

alpha <- 0.05   # significance level

#m=29  # for real data of "phenotype_reduced.csv"

# sample variance of a measure vs number of samples
val <- NULL  
##########for(n in seq(10, 500, 10)) {   
for(n in seq(20, 500, 50)) {   
##################for(n in seq(300, 300, 10)) {   # to draw histogram of density
   # perform 'p' trials to evaluate sample variance
   ans <- NULL
   for(i in 1:p) {
      x <- matrix(, n, m) 
      #print("OK1")

      # select random or correlated case...
      if(oracle==0) {
         # generate random phenotypes 'x' of n samples and m phenotypes
         for(j in 1:m) { x[,j] <- rnorm(n) }
      } else {
         # generate correlated phenotypes 'x' of n samples and m phenotypes
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

      # test of correlation with a matrix
      # rcorr() computes a matrix of Pearson's r (or Spearman's rho rank 
      # correlation coefficients) for all possible pairs of columns of a matrix
      # It also returns p-values.
      #re <- rcorr(t(npheno), type="spearman")  
      #re <- rcorr(t(npheno), type="pearson")  
      #p_vec <- re$P[lower.tri(re$P)]   # vector of p-values re$P

      # test of correlation using cor.test() for Kendall's statistics
      p_vec <- NULL
      for(a in 1:(n-1)) {
         for(b in (a+1):n) {
            pp <- cor.test(npheno[a,], npheno[b,], method="pearson")
      	    #pp <- cor.test(npheno[a,], npheno[b,], method="spearman")
      	    #pp <- cor.test(npheno[a,], npheno[b,], method="kendall")
      	    # print(pp$p.value)
      	    p_vec <- c(p_vec, pp$p.value)
         }
      }
      # end of test

      res <- cgraph(p_vec, n, alpha)  # construct a graph and graph density
      ans <- c(ans, res[[1]])    # graph density
      #ans <- c(ans, res[[2]])    # graph transitivity (or clustering coeff)
   }
   #val <- rbind(val, c(n, var(ans)))                    # old: variance only
   val <- rbind(val, c(n, mean(ans), sd(ans), var(ans)))   # new: mean and SD
   ########return(ans)    # to draw histogram of density
}

if(oracle==0) {
   write.csv(val, "gden_random_pearson.csv")
   #write.csv(val, "gden_random_spearman.csv")
   #write.csv(val, "out3_random.csv")
   #write.csv(val, "out3_random_trans.csv")
   #write.csv(val, "out3_random_pearson.csv")
   #write.csv(val, "out3_random_kendall.csv")
   #write.csv(ans, "density_random.csv")  # to draw histogram of density
} else {
   write.csv(val, "gden_corr_pearson.csv")
   #write.csv(val, "gden_corr_spearman.csv")
   #write.csv(val, "out3_corr.csv")
   #write.csv(val, "out3_corr_trans.csv")
   #write.csv(val, "out3_corr_kendall.csv")
   #write.csv(ans, "density_corr.csv")  # to draw histogram of density
}

#################################################
#Rscript Density3.R  30  20  0 &    # for random
#Rscript Density3.R  30  20  1 &    # for correlated
#################################################################

# to find frequency....  11/25/2021
#Rscript Density3.R  30  100  0 &    # for random
#Rscript Density3.R  30  100  1 &    # for correlated

# to draw histogram of density with 200 samples
#Rscript Density3.R  30  200  0 &    # for random
#Rscript Density3.R  30  200  1 &    # for correlated
