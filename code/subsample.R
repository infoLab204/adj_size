# select n subsamples from original max nn samples

subsample <- function(geno, pheno1, pheno2, s, m) {
  # s is number of subsamples
  # m is number of phenotypes

  nn <- length(geno)                  # nn is max sample size
  idx <- sample(1:nn, s, replace=F)   # s is subsample size 
  idx <- sort(idx)

  # subsampling s samples
  sgeno <- NULL
  spheno1 <- matrix(, s, m) 
  spheno2 <- matrix(, s, m) 

  for(i in 1:s) {
    sgeno <- c(sgeno, geno[idx[i]])
    spheno1[i, ] <- pheno1[idx[i], ]   # independent phenotype
    spheno2[i, ] <- pheno2[idx[i], ]   # dependent phenotype
  }

  # return subsamples 
  return(list(sgeno, spheno1, spheno2))   
}
