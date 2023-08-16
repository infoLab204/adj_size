// first approximation for lambda

double flambda(int *x, double *y, int n, double lambda, double alpha) 
{
  double sum, fi, df, pp ;
  int i ;

  for(sum=0.0, i=0 ; i < n ; i++) {
    pp=pow(x[i], 1.0-alpha) ;
    fi=lambda/(1.0-alpha)*(1.0-pp)+log(lambda)-alpha*log(x[i]) ;
    df=1.0/lambda+(1.0-pp)/(1.0-alpha) ; 
    sum+=(log(y[i])-fi)*df ;
  }
  return sum ;
}

  
