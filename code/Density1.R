#################################################################
# excutable R-script for evaluating sample variance of graph density in 
# real and shuffled data of 44K, including Yesan NAM data.
# Note: shffled data can be obtained using Shuffle.R in current directory 
#################################################################
# include packages
library(Hmisc)   # rcorr() 
library(igraph)
library(MASS)
library(ppcor)   # for partial correlation coefficient

# estimate partial correlation coefficient and graphical analysis
source("eval_corr.R") # estimate correlation coefficient
source("eval_pcorr.R") # estimate partial correlation 
source("shuffle.R")   # shuffling data
source("cgraph.R")   # graphical analysis of data

args <- commandArgs(TRUE)

p <- as.numeric(args[1])  # number of trials to estimate sample variance
oracle <- as.numeric(args[2]) # select either shuffled or original data

if(oracle==0) { ########## shuffled data from Shuffle.R
   #pheno <- read.csv("phenotype_reduced_shuffled.csv", header=T) 
   #pheno <- read.csv("phenotype_30traits_shuffled.csv", header=T) 
   pheno <- read.csv("phenotype_nam_7traits_shuffled.csv", header=T) 
   pheno <- pheno[, -c(1,2)]  # exclude the first column of sample ID
} else{  ############# real data 
   #pheno <- read.csv("phenotype_reduced.csv", header=T) 
   #pheno <- read.csv("phenotype_30traits.csv", header=T) 
   pheno <- read.csv("phenotype_nam_7traits.csv", header=T) # Yesan NAM data 
   pheno <- pheno[, -1]  # exclude the first column of sample ID
}

alpha=0.05   # significance level

x <- as.matrix(pheno)   # check matrix type of pheno
nn=nrow(x)   # number of samples
m=ncol(x)   # number of phenotypes

# iteration for variance vs number of samples 
val <- NULL
for(n in seq(20, nn, 50)) {
######for(n in seq(10, nn, 10)) {
#####for(n in seq(150, nn, 10)) {   ############# to draw histogram of density
   ans <- NULL
   # perform 'p' iterations to evaluate sample variance
   for(i in 1:p) {
      # sampling 
      idx <- sample(1:nn, n, replace=F)
      spheno <- x[idx, ]

      # normalize pheno "x" before evaluationg p-values 
      # npheno is a matrix class
      npheno <- NULL   # normalized phenotypes
      for(j in 1:m) {
         y <- (spheno[, j] - mean(spheno[, j]))/sd(spheno[, j])
      	 npheno <- cbind(npheno, y)
      }

      ################# test of correlation with a matrix
      # test for correlation using 'Hmisc' package
      # rcorr() computes a matrix of Pearson's r (or Spearman's rho rank 
      # correlation coefficients) for all possible pairs of columns of a matrix
      # It also returns p-values.
      #re <- rcorr(t(npheno), type="spearman")  
      #re <- rcorr(t(npheno), type="pearson")  
      #p_vec <- re$P[lower.tri(re$P)]   # vector of p-values re$P

      ###### test of correlation using cor.test() for Kendall's statistics
      p_vec <- NULL
      for(a in 1:(n-1)) {
         for(b in (a+1):n) {
      	    #pp <- cor.test(npheno[a,], npheno[b,], method="spearman")
      	    pp <- cor.test(npheno[a,], npheno[b,], method="kendall")
            #pp <- cor.test(npheno[a,], npheno[b,], method="pearson")
      	    # print(pp$p.value)
      	    p_vec <- c(p_vec, pp$p.value)
        }
      }
      ################### end of test

      re <- cgraph(p_vec, n, alpha)  # construct graph and its density
      ans <- c(ans, re[[1]])    # graph density
      #ans <- c(ans, re[[2]])     # graph transitivity
   }
   #val <- rbind(val, c(n, var(ans)))   # old: variance only
   val <- rbind(val, c(n, mean(ans), sd(ans), var(ans)))  # new: mean and sd
   ####return(ans)  ###################  to draw histgram of density
}
#var(ans)

if(oracle==0){
   write.csv(val, "gden_nam_shuffled_kendall.csv")
   #write.csv(val, "gden_nam_shuffled_spearman.csv")
   #write.csv(val, "gden_rshuffled_pearson.csv")
   #write.csv(val, "gden_rshuffled_kendall.csv")
   #write.csv(val, "gden_rshuffled_spearman.csv")
   #write.csv(val, "out1_shuffled_spearman.csv")
   #write.csv(val, "out1_shuffled_pearson.csv")
   #write.csv(val, "out1_shuffled_kendall.csv")
   #write.csv(ans, "density_real_shuffle.csv")  # to draw histogram of density
} else{
   write.csv(val, "gden_nam_kendall.csv")
   #write.csv(val, "gden_nam_spearman.csv")
   #write.csv(val, "gden_real_pearson.csv")
   #write.csv(val, "gden_real_kendall.csv")
   #write.csv(val, "gden_real_spearman.csv")
   #write.csv(val, "out1_real_spearman.csv")
   #write.csv(val, "out1_real_trans.csv")
   #write.csv(val, "out1_real_pearson.csv")
   #write.csv(val, "out1_real_kendal.csv")
   #write.csv(ans, "density_real_unshuffle.csv")  # to draw histogram of density
}

#################################################################
#Rscript Density1.R 20 0 &     # for shuffle
#Rscript Density1.R 20 1 &     # for unshuffle
#################################################################


######################### for graph density ###########################
#Rscript Density1.R 200 0 &     # for shuffle
#Rscript Density1.R 200 1 &     # for unshuffle
#################################################################
