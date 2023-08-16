#################################################################
# excutable R-script for evaluating sample variance of graph density in 
# real and shuffled data of 580K categorical data
# Note: shffled data can be obtained using Shuffle.R in current directory 
#################################################################
# include packages
#library(Hmisc)   # rcorr() 
library(igraph)
library(MASS)
#library(ppcor)   # for partial correlation coefficient

# estimate partial correlation coefficient and graphical analysis
#source("eval_corr.R") # estimate correlation coefficient
#source("eval_pcorr.R") # estimate partial correlation 
#source("shuffle.R")   # shuffling data
source("cgraph.R")   # graphical analysis of data

args <- commandArgs(TRUE)

p <- as.numeric(args[1])  # number of trials to estimate variance of density
oracle <- as.numeric(args[2]) # select either shuffled or original data

if(oracle==0) { ########## shuffled data using Shuffle.R
   #pheno <- read.csv("phenotype_30traits_shuffled.csv", header=T) 
   pheno <- read.csv("phenotype_580k_14qtraits_s.csv", header=F) 
   pheno <- pheno[, -c(1,2)]  # exclude the first column of sample ID
} else{  ############# real data 
   #pheno <- read.csv("phenotype_30traits.csv", header=T) 
   pheno <- read.csv("phenotype_580k_14qtraits.csv", header=F) # Categorical 
   pheno <- pheno[, -1]  # exclude the first column of sample ID
}

alpha=0.05   # significance level

x <- as.matrix(pheno)   # check matrix type of pheno

nn=nrow(x)   # number of samples (accession)
m=ncol(x)   # number of categorical phenotypes

# iteration to estimate variance vs number of samples 
val <- NULL
for(n in seq(10, nn, 30)) {
#for(n in seq(10, 50, 20)) {
   ans <- NULL
   # perform 'p' iterations to evaluate variance

   for(i in 1:p) {
      # sampling 
      idx <- sample(1:nn, n, replace=F)
      npheno <- x[idx, ]

      ###### test of correlation using non-parametric sign test
      p_vec <- NULL
      for(a in 1:(n-1)) {
         for(b in (a+1):n) {
      	    #pp <- cor.test(npheno[a,], npheno[b,], method="kendall")
	    #pp <- chisq.test(npheno[a,], npheno[b,])

	    # count non-equal categorical phenotypes out of m
	    neq <- 0
	    for(k in 1:m) { if(npheno[a,k]==npheno[b,k]) neq <- neq+1 }
	    pp <- binom.test(neq, m)
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

if(oracle==0){
   write.csv(val, "gden_580k_sign_shuffled.csv")
   #write.csv(val, "gden_nam_shuffled_spearman.csv")
   #write.csv(ans, "density_real_shuffle.csv")  # to draw histogram of density
} else{
   write.csv(val, "gden_580k_sign.csv")
   #write.csv(val, "gden_nam_spearman.csv")
   #write.csv(ans, "density_real_unshuffle.csv")  # to draw histogram of density
}

#################################################################
#Rscript Density2.R 20 0 &     # for shuffled
#Rscript Density2.R 20 1 &     # for unshuffled
#################################################################

