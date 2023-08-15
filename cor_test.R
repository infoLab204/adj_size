# estimate the variance of graph density
cor_test <- function(input, p) {
  # input contains genotype vector and phenotype matrix
  geno <- input[[1]]    # genotype vector of size nn
  pheno1 <- input[[2]]  # nn*m independent phenotype matrix 
  pheno2 <- input[[3]]  # nn*m dependent phenotype data
  nn <- length(geno)   # number of samples 
  m <- ncol(pheno1)    # number of phenotypes
  # p is number of indepedent trials

  alpha <- 0.05   # significance level
  val1 <- NULL  # sample variance of a measure vs number of samples
  val2 <- NULL  # sample variance of a measure vs number of samples
  ###### max subsample size is 500, which is less than nn=1000
  ###### 500 < nn, because we like to maintain about randomness of sampling
  ###### for(n in seq(10, 500, 10)) {   # vary subsample size n, original
  #for(n in seq(20, 500, 50)) {   # vary subsample size n, plot density
  for(n in seq(30, 500, 30)) {   # vary subsample size n, plot density
    ans1 <- NULL ; ans2 <- NULL
    for(i in 1:p) {
      # select n samples out of nn
      sub <- subsample(geno, pheno1, pheno2, n, m)
      sgeno <- sub[[1]]    # subsampled genotype
      spheno1 <- sub[[2]]  # subsampled independent phenotype
      spheno2 <- sub[[3]]  # subsampled dependent phenotype

      ####### In the following, we have to select correlation test statistis
      ################# select correlation test: Pearson, Spearman, Kendall
      # test of correlation with a matrix npheno
      # rcorr() computes a matrix of Pearson's r (or Spearman's rho rank 
      # correlation coefficients) for all possible pairs of columns of a matrix
      # It also returns p-values.

      #re <- rcorr(t(spheno1), type="spearman")  #### for spearman's test
      #re <- rcorr(t(spheno1), type="pearson")    #### for pearson's test
      #p_vec <- re$P[lower.tri(re$P)]   # vector of p-values re$P

      ########### test of correlation using cor.test() for Kendall's statistics
      p_vec1 <- NULL ;       p_vec2 <- NULL
      for(a in 1:(n-1)) {
         for(b in (a+1):n) {
           #pp1 <- cor.test(spheno1[a,], spheno1[b,], method="spearman") # indep
           #pp2 <- cor.test(spheno2[a,], spheno2[b,], method="spearman") # dep
           pp1 <- cor.test(spheno1[a,], spheno1[b,], method="pearson") # indep
           pp2 <- cor.test(spheno2[a,], spheno2[b,], method="pearson") # dep
           #pp1 <- cor.test(spheno1[a,], spheno1[b,], method="kendall")
           #pp2 <- cor.test(spheno2[a,], spheno2[b,], method="kendall")
           #print(pp1$p.value)
      	   p_vec1 <- c(p_vec1, pp1$p.value)
      	   p_vec2 <- c(p_vec2, pp2$p.value)
         }
      }
      ############# end of (Kendall) test

      res1 <- cgraph(p_vec1, n, alpha)  # construct a graph and graph density
      ans1 <- c(ans1, res1[[1]])     # graph density
      #ans1 <- c(ans1, res1[[2]])    # graph transitivity (or clustering coeff)
 
      res2 <- cgraph(p_vec2, n, alpha)  # construct a graph and graph density
      ans2 <- c(ans2, res2[[1]])     # graph density
      #ans2 <- c(ans2, res2[[2]])    # graph transitivity (or clustering coeff)
    }
    val1 <- rbind(val1, c(n, mean(ans1), sd(ans1), var(ans1)))
    val2 <- rbind(val2, c(n, mean(ans2), sd(ans2), var(ans2)))
    #val1 <- rbind(val1, c(n, var(ans1)))
    #val2 <- rbind(val2, c(n, var(ans2)))
    #######return(ans1)    #################### to draw histogram of density
  }

  ####### for both variance of graph density and the mean, sd of graph density
  write.csv(val1, "gden_random_pearson.csv")
  write.csv(val2, "gden_corr_pearson.csv")
  #write.csv(val1, "gden_random_kendall.csv")
  #write.csv(val2, "gden_corr_kendall.csv")
  #write.csv(val1, "gden_random_spearman.csv")
  #write.csv(val2, "gden_corr_spearman.csv")

  #write.csv(val1, "density_indep_spearman.csv")   # variance of graph density
  #write.csv(val2, "density_dep_spearman.csv")
  #write.csv(val1, "density_indep_pearson.csv")
  #write.csv(val2, "density_dep_pearson.csv")
  #write.csv(val1, "density_indep_kendall.csv")
  #write.csv(val2, "density_dep_kendall.csv")


  #write.csv(ans1, "histo_indep_spearman.csv") #  histogram of density
  #write.csv(ans2, "histo_dep_spearman.csv") #  histogram of density

  ################## old #############################
  #write.csv(val, "density_random_kendall.csv")
  #write.csv(ans, "histo_random_pearson.csv") # to draw histogram of density
  #write.csv(val, "density_corr_pearson.csv")
  #write.csv(val, "density_corr_spearman.csv")
  #write.csv(val, "density_corr_kendall.csv") 
  #write.csv(ans, "histo_corr_pearson.csv") # to draw histogram of density
  #write.csv(ans, "histo_corr_spearman_test.csv") # to draw histogram of density
}