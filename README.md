# Adj_size: a computational model for estimating adjusted sample size     
Sunhee Kim and Chang-Yong Lee    

Adj_size represents R and C scripts to perform various hypothesis tests and estimate the parameters of the graph density model to estimate the adjusted sample size using simulated and real data sets.    

We proposed a computational model and method to account for sample dependence in hypothesis testing. We used graph theory to quantify the degree of dependence in terms of graph density and model the variance of graph density as a function of sample size. We simulated model-free and genotype-phenotype association data sets to support the validity of the proposed method. We then applied the method to rice genomic data sets to estimate the adjusted sample size as a guide to account for sample dependence in hypothesis testing. We provide the R and C scripts together with real data sets in order for the readers to reproduce the results discussed in the manuscript.    

## Data sets
*	preprocessed 44K_SNP data
*	preprocessed 580K_SNP data



## R script for estimating p values 
### 1. Generating model-free data and estimating p values: Pvalue_nomodel.R       
*	eval_corr.R    
*	usage: Rscript Pvalue_nomodel.R  20  20  0 &    # for independent data    
*	usage: Rscript Pvalue_nomodel.R  20  20  1 &    # for dependent data    

### 2. Generating model-based data and estimating p values: Pvalue_model.R     
*	gene_input.R: generate input    
*	subsample.R: subsampling     
*	tf_test.R: t- and F-test   
*	usage: Rscript Pvalue_model.R  1100  20  0.3  0.4  &     



## R script for graph density modeling
### 1.	Estimating variance of graph density: model-free data      
●	Density3.R: main script    
●	eval_corr.R: estimate correlation coefficient    
●	cgraph.R: graphical analysis of data    

### 2.	Estimating variance of graph density: model-based data    
●	Gdensity.R: main script    
●	gene_input.R: generate input    
●	subsample.R: subsampling   
●	cgraph.R: graphical analysis of data    
●	cor_test.R: estimate correlation coefficient    

### 3.	Estimating variance of graph density: 44K_SNP data    
●	Density1.R: main script    
●	eval_corr.R: estimate correlation coefficient    
●	cgraph.R: graphical analysis of data    
●	shuffle.R: shuffling data    

### 4.	Estimating variance of graph density: 580K_SNP data    
●	Density2.R: main script    
●	cgraph.R: graphical analysis of data   

## C codes for modeling variance of graph density and estimating parameters

### 1.	Estimate parameters: least_square.c        
*	flambda.c: estimate parameter lambda    
*	falpha.c: estimate parameter alpha    

### 2.	Data fitting with parameters: fit.c    

### 3.	Estimate adjusted sample size: solveg.c    

