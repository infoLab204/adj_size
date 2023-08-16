# t-test and F-test for independent and dependent data
# varies sample size while fixing max sample, correlation, and phenotypes

tf_test <- function(input) { 
  # input contains genotype vector and phenotype matrix
  geno <- input[[1]]    # genotype vector of size nn
  pheno1 <- input[[2]]  # nn*m independent phenotype matrix 
  pheno2 <- input[[3]]  # nn*m dependent phenotype data

  nn <- length(geno)   # max number of genotypes
  m <- ncol(pheno1)    # number of phenotypes

  # save test results over m phenotypes
  # ncol=2 means p-values for two cases (independent and dependent data)
  htest <- matrix(data=NA, nrow=m, ncol=2)   

  # save 'p' independent test results averaged over m phenotypes
  p <- 20  # p is number of independent trials        
  
  # the first row of 'result' is reserved for the subsample size
  # value 'ncol' should be int(nn/100)*2  (eg) int(1100/100)*2=20
  # model1
  result <- matrix(data=NA, nrow=p+1, ncol=as.integer(nn/100)*2) 

  # model2: NOT used
  #sn <- 510  # model2, max subsample size
  #result <- matrix(data=NA, nrow=p+1, ncol=as.integer(sn/50)*2) 

  # for each subsample size 'n' with index 'id'
  id=1 

  # select model
  for(n in seq(from=100, to=nn, by=100)) {   # model1 vary subsample size n 
  #for(n in seq(50, sn, 50)) {   # model2, vary subsample size n. NOT used
    # enter subsample size 'n' into the first row
    result[1, 2*id-1] <- n
    result[1, 2*id] <- n

    # p independent trials to obtain statistics
    for(j in 1:p) {
      # subsampling of size n, (eg: 100 <= n <= 1100)
      # we fixed max sample size (eg, 1100) and sampling n samples out of max
      # The reason is that the random correlation matrix fluctuates...

      # select n samples out of nn
      sub <- subsample(geno, pheno1, pheno2, n, m)

      sgeno <- sub[[1]]    # subsampled genotype
      spheno1 <- sub[[2]]  # subsampled independent phenotype
      spheno2 <- sub[[3]]  # subsampled dependent phenotype

      # perform t- or F-test given each of m phenotypes
      for(i in 1:m) {
	####### we should make a choice: either F- or t-test 
	#######  F-test
	fpheno1 <- spheno1[,i]
	fpheno2 <- spheno2[,i]
	df1 <- data.frame(sgeno, fpheno1)
	df2 <- data.frame(sgeno, fpheno2)

	res1 <- summary(aov(fpheno1 ~ sgeno, df1))
	res2 <- summary(aov(fpheno2 ~ sgeno, df2))

	htest[i, 1] <- res1[[1]][[5]][1]  # p-value from independent data
	htest[i, 2] <- res2[[1]][[5]][1]  # p-value from dependent data

      	######## t-test
      	#htest[i, 1] <- cor.test(sgeno, spheno1[, i])$p.value
      	#htest[i, 2] <- cor.test(sgeno, spheno2[, i])$p.value
	###### end of test
      }

      # htest contains two p-values of independent and dependent data 
      # mean p-value over m phenotypes
      ans <- apply(htest, 2, mean)

      result[j+1, 2*id-1] <- ans[1]   # independent data
      result[j+1, 2*id] <- ans[2]     # dependent data
    } # end of p independent trials

    id <- id+1   
  } # end of each subsample 'n'

  #write.csv(result, "pvalue_t_test.csv")  ###  contain both indep and dep
  #write.csv(result, "pvalue_F_test.csv")
  write.csv(result, "pvalue_F_test.csv")

  #return(result)
}
