# generate input data: genotype vector and shuffled phenotype matrix
# according to the model of genotype-phenotype

gene_input <- function(n, m, q, kappa)
{
   # n is max number of samples, m is number of phenotypes for each genotype
   # q is maf, kappa is shuffling rate

   # generate m uniform random numbers for m phenotypes 
   pi <- runif(m)   
   # generate normal random vector for the evaluation of each phenotype
   epsilon <- rnorm(n)   # epsilon ~ N(0, 1)

   # generate n genotypes with maf q=0.3
   geno <- sample(0:2, n, replace=TRUE, prob=c((1-q)*(1-q), 2*q*(1-q), q*q))

   # generate m independent phenotypes according to the model
   z <- matrix(, n, m)   # n*m independent phenotype matrix 
   for(j in 1:m) 
     z[,j]<-sqrt(1-pi[j])*epsilon+sqrt(pi[j]/(2*q*(1-q)))*geno

   # shuffle phenotypes to create non-associated genotype-phenotype pairs
   sub <- as.integer(kappa * n)    # number of shuffled phenotypes
   sub <- ifelse(sub%%2==0, sub, sub+1)  # making even number for swap

   # swapping for partial association
   for(j in 1:m) {
     idx <- sample(1:n, sub, replace=F) 
     for(k in seq(1, sub, 2)) swap(z[idx[k], j], z[idx[k+1], j])
   }
   ### end of generating independent phenotype data

   # generate n*n correlation matrix, where n is the max number of samples 
   corr <- rcorrmatrix(n, alphad=1)   # random correlation matrix
   # making positive definite & Cholesky decompo.
   nchol <- chol(nearPD(corr)) 

   #write.table("inside input1", "test1.txt")

   # generate dependent phenotype data
   x <- matrix(, n, m)   # n*m dependent phenotype matrix 
   for(j in 1:m) 
     x[,j]<-sqrt(1-pi[j])*nchol%*%z[,j]+sqrt(pi[j]/(2*q*(1-q)))*geno

   #write.table("inside input2", "test1.txt", append=TRUE)

   # return genotype, indpendent, and dependent phenotype data
   return(list(geno, z, x))
}
 