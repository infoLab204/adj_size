# shuffle given data to estimate correlation coefficient of random pheno
# shuffling a data.frame by rows of each column
# returns a shuffled data
shuffle <- function(d) {
  df <- NULL
  n <- nrow(d)
  m <- ncol(d)
  for(i in 1:m) {
    s <- d[sample(n, n), i]
    df <- cbind(df, s) 
  }
  return(as.data.frame(df))
}