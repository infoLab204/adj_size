# excutable R-script for the model of association between genotype and phenotype
# generate independent and dependent data to perform t- and F-test
# plot p-values against sample size

# include packages
library(MASS) 
library(seqinr)   # for swap()
library(clusterGeneration)   # for rcorrmatrix(), random correlation matrix 
library(lmf)   	# package "lmf" for nearPD() for covariance matrix

# generate input: correlation coefficient, genotypes, and random interaction 
source("gene_input.R") # generate input
source("subsample.R")  # subsampling 
source("tf_test.R")    # F- and t-test 

args <- commandArgs(TRUE)

# args[1] is number of genotypes
# args[2] is phenotype size
# args[3] is maf
# args[4] is kappa, shuffling rate
# generate input data according to the model
input <- gene_input(as.numeric(args[1]), as.double(args[2]), as.numeric(args[3]), as.double(args[4]))

# association test using input data
res <- tf_test(input)

# save results
# write.csv(res, "pvalue_F_test.csv")  
# write.csv(res, "pvalue_t_test.csv")

####################### command #################################
#Rscript Pvalue_model.R  1100  20  0.3  0.4  &   
#################################################################
