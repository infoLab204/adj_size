// fit the model with estimated parameters
// a.out >  output3_random_pearson
// a.out > output1_real_spearman

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv)
{
  int n;
  double pp, fi, alpha, lambda ;

  // estimated parameters of simulated data
  alpha=1.987619 ;    // simulated data, random spearman
  lambda=3.844194 ; 

  alpha=2.066264;     // simulated data, correlated spearman
  lambda=1.704572 ;

  alpha=1.992540 ;    //  simulated data, random Kendall
  lambda=3.776760 ;   

  alpha=2.106959 ;    //  simulated data, corr kendall
  lambda=1.118714 ;   


  alpha=2.093196 ;  // corr spearman scale=5.0
  lambda=1.072395 ; // corr spearman scale=5.0

  alpha=2.069634 ;  // corr spearman scale=3.0
  lambda=1.185966 ; // corr spearman scale=3.0

  alpha=2.072132 ;  // corr spearman scale=10.0
  lambda=1.160616 ; // corr spearman scale=10.0

  alpha=2.301045 ;  // corr spearman scale=20.0
  lambda=1.160616 ; // corr spearman scale=20.0

  alpha=2.106756 ;  // corr spearman scale=7.0
  lambda=1.227719  ; // corr spearman scale=7.0

  /// --------------------------------------------------------
  alpha=1.490526; // corr1_trans
  lambda=0.694895; 

  // alpha=1.605454; // random1_trans
  //lambda=0.757132; 
  ///--------------- real data -------------------------//
  // simulated data of same size as real data, random spearman, not used 
  alpha=2.153;   
  lambda=3.318;

  // simulated data of same size as real data, correlated spearman, not used
  alpha=1.055078;  
  lambda=1.486663; 


  alpha=1.682557 ;    // real kendall 
  lambda=0.822935  ;  // real kendall

  alpha=2.506600 ;    // real shuffled kendall
  lambda=0.828457 ;   // real shuffled kendall

  alpha=1.666006 ;    // real random kendall 
  lambda=3.709132  ;  // real random kendall

  alpha=1.547746 ;    // real spearman
  lambda=0.840864 ;   // real spearman

  alpha=2.526121 ;    // real shuffled spearman
  lambda=1.499796  ;  // real shuffled spearman

  alpha=2.152928 ;    // real random spearman
  lambda=3.318316 ;   // real random spearman

  alpha=1.690494 ;    // real pearson
  lambda=0.821952  ;  // real pearson

  alpha=1.786100 ;    // real random pearson
  lambda=3.844611  ;  // real random pearson

  alpha=2.569141 ;    // real shuffle pearson
  lambda=1.557751  ;  // real shuffle pearson
  /////////////////////////////////////////////////////

  alpha=1.526899 ;  // real spearman
  lambda=0.760727 ;

  alpha=2.232351 ;  // shuffled spearman
  lambda=3.558407 ; 

  alpha=1.076775 ;  // shuffled pearson
  lambda=1.954146 ;

  alpha=1.577892 ;  // real pearson
  lambda=0.746412 ;

  alpha=1.605463 ;  // real kendall
  lambda=0.738669 ;

  /*
  alpha=2.186375 ;  // shuffled kendall
  lambda=3.866892 ;
  */


  alpha=1.990261 ;   // nomodel random pearson
  lambda=3.924713 ; 


  alpha=1.523790 ;   // nomodel corr pearson
  lambda=2.799959  ;

  alpha=0.524693 ;  // 580k real 14q
  lambda=0.091047 ; 


  alpha=1.662493 ;  // 580k real 14q shuffled
  lambda=2.209005 ;

  //for(n=1 ; n < 235 ;  n++) {   // for real data (44K SNP)
  //for(n=1 ; n < 550 ;  n++) {   // for simulated data
  for(n=1 ; n < 850 ;  n+=2) {   // for 850k 
    pp=pow(1.0*n, 1.0-alpha) ;
    fi=lambda/(1.0-alpha)*(1.0-pp)+log(lambda)-alpha*log(1.0*n) ;
    printf("%d\t%.10lf\n", n, exp(fi)) ;
  }
  return 0 ;
}

