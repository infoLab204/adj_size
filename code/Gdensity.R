#################################################################
# excutable R-script for evaluating sample variance of graph density and 
# transitivity in correlated and uncorrelated values (or phenotypes) by using
# cor.test instead of rcorr for Kendall statistics
#################################################################
# include packages
library(Hmisc)   # rcorr(), correlation test
library(igraph)
library(MASS)
library(seqinr)  # for swap()
library(clusterGeneration)   # for rcorrmatrix() 
library(lmf)   	# package "lmf" for nearPD() for covariance matrix

# generate input: correlation coefficient, genotypes, and random interaction 
source("gene_input.R") # generate input
source("subsample.R")  # subsampling 
source("cgraph.R")   # graphical analysis of data
source("cor_test.R") # estimate correlation coefficient

args <- commandArgs(TRUE)

# args[1] is max sample size (eg) 1000
# args[2] is number of phenotypes  (eg) 30
# args[3] is maf  (eg) 0.3
# args[4] is kappa, shuffling rate  (eg) 0.4
# args[5] is number of trials to estimate sample variance (eg) 20

# generate input data according to the model
input <- gene_input(as.numeric(args[1]), as.double(args[2]), as.numeric(args[3]), as.double(args[4]))

# correlation test to construct a graph and estimate graph density
cor_test(input, as.numeric(args[5])) 

################################################
# for the variance of graph density
#Rscript Gdensity.R  1000 30  0.3  0.4  20 &    # should be 1000
#################################################################

################################################
# for histogram of density, use 300 samples and 200 independent trials
#Rscript Gdensity.R  300 30  0.3  0.4  200 &    
#################################################################


