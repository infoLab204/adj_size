// first approximation for alpha

double falpha(int *x, double *y, int n, double alpha, double lambda)
{
  double sum, pp, fi, df ;
  int i ;

  sum=0.0 ; 
  for(i=0 ; i < n ; i++) {
    pp=pow(x[i], 1.0-alpha) ;
    fi=lambda/(1.0-alpha)*(1.0-pp)+log(lambda)-alpha*log(x[i]) ;
    df=lambda/((1.0-alpha)*(1.0-alpha))*(1.0-pp)+(lambda*pp/(1.0-alpha)-1.0)*log(x[i]) ;
    sum+=(log(y[i])-fi)*df ;
  }
  return sum ; 
}
