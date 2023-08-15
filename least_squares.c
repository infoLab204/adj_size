/*                                                                            
 * Least squares estimation by using bisection method for solving 
 *      simulateous equations for lambda and alpha      
 *         
 * Estimate alpha and lambda for upper and lower bounds, and take average
 *
 * $a.out  input1_real_spearman    outpara1_real_spearman   1.0e-5  0.001
 * $a.out  input1_rshuffled_spearman   outpara1_rshuffled_spearman   1.0e-5  0.001
 *
 * $a.out  input_random_kendall   outpara_random_kendall   1.0e-5  0.0001
 *
 * $a.out  input_random1_trans  outpara_random1_trans  1.0e-5  0.0001
 * $a.out  input_corr1_trans  outpara_corr1_trans  1.0e-5  0.0001
 *
 *   input: input_random, or input_corr, or gden_real

 * $a.out  input3_random_pearson   outpara3_random_pearson   1.0e-5  0.0001
 * $a.out  input3_corr_pearson   outpara3_corr_pearson   1.0e-5  0.0001
 *
 * $a.out input2_580k outpara2_580k 1.0e-05 0.0001
 *  
 */

#define JMAX 40   // maximum allowed number of bisections  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double flambda(int *x, double *y, int n, double lambda, double alpha)  ;
double falpha(int *x, double *y, int n, double alpha, double lambda) ;
double rtbis(double (*func)(int*, double*, int, double, double), double x1, double x2, double para, double xacc, int *dax, double *day, int n) ;

#include "flambda.c"
#include "falpha.c"

int main(int argc, char **argv)
{
  FILE *fi, *fs, *fo ;
  int *dx, n, i, val, chr, count ; 
  double *dy, a1, a2, b1, b2, alpha, lambda, oalpha, olambda, acc, valy ;
  double (*pfun1)(int *x, double *y, int n, double a, double l) ;
  double (*pfun2)(int *x, double *y, int n, double a, double l) ;

  if(argc!=5) {
    printf("enter correct parameters...\n") ;
    exit(1) ;
  }

  fi=fopen(argv[1], "r") ;   // input data 
  fo=fopen(argv[2], "w") ;
  acc=atof(argv[3])  ;  // error or accuracy of the bisection method   
  alpha=atof(argv[4]) ;        

  pfun1=flambda ;
  pfun2=falpha ;

  //---------------------------------------------
  //a1=0.0000001 ;  // lower bound of lambda for input_corr1_trans
  //b1=10.0 ;     // upper bound of lambda for input_corr1_trans

  //a1=0.0000001 ;  // lower bound of lambda for input_random1_trans
  //b1=10.0 ;     // upper bound of lambda for input_random1_trans

  //a1=0.0000001 ;  // lower bound of lambda for input_random1  (or kendall)
  //b1=10.0 ;     // upper bound of lambda for input_random1  (or kendall)

  //a1=0.00001 ;  // lower bound of lambda for input_corr1
  //b1=100.0 ;    // upper bound of lambda for input_corr1

  a1=0.0000001 ;  // lower bound of lambda for input_corr_kendall
  b1=100.0 ;     // upper bound of lambda for input_corr_kendall and correlation with scale=5.0, scale=3.0, scale=7.0, scale=10.0  spearman

  a1=0.00000001 ;  // lower bound of lambda for input_corr_kendall
  b1=1000.0 ;     // upper bound of lambda for input_corr_kendall and correlation with scale=20 spearman

  //-----------------------for real data -------------------
  a1=0.00001 ;  // lower bound of lambda for input_random and kendall
  b1=10.0 ;     // upper bound of lambda for input_random

  a1=0.00001 ;  // lower bound of lambda for input_pearson
  b1=10.0 ;     // upper bound of lambda for input_pearson

  //a1=0.01 ;  // lower bound of lambda for input_corr or input_real
  //b1=1.0 ;   // upper bound of lambda for input_corr or input_real

  // -------------------------------------------------
  // -------------------------------------------------

  //a2=0.0000001 ;  // lower bound of alpha for input_corr1_trans
  //b2=10.0 ;     // upper bound of alpha for input_corr1_trans

  //a2=0.0000001 ;  // lower bound of alpha for input_random1_trans
  //b2=10.0 ;     // upper bound of alpha for input_random1_trans

  //a2=0.0000001 ;  // lower bound of alpha for input_random1 (or kendall)
  //b2=10.0 ;     // upper bound of alpha for input_random1  (or kendall)

  //a2=0.00001 ;  // lower bound of alpha for input_corr1
  //b2=100.0 ;     // upper bound of alpha for input_corr1

  a2=0.0000001;  // lower bound of alpha for input_corr_kendall
  b2=100.0;     // upper bound of alpha for input_corr_kendall and scale=5.0, scale=3.0, scale=7, scale=10.0  spearman

  a2=0.00000001;  // lower bound of alpha for input_corr_kendall
  b2=1000.0;     // upper bound of alpha for input_corr_kendall and scale=5.0, scale=20  spearman

  //--------------------------- for real data -----------------------
  a2=0.00001 ;  // lower bound of alpha for input_random and kendall
  b2=10.0 ;     // upper bound of alpha for input_random

  a2=0.001 ;  // lower bound of lambda for input_pearson
  b2=10.0 ;     // upper bound of lambda for input_pearson

  //a2=0.01 ;  // lower bound of alpha for input_corr or input_real
  //b2=1.0 ;   // upper bound of alpha for input_corr or input_real


  a1=0.0000001 ;  // lower bound of lambda for input_real_shu (44K SNP)
  b1=10.0 ;     // upper bound of lambda for input_real_shu

  a2=0.0000001 ;  // lower bound of alpha for input_real_shu
  b2=10.0 ;     // upper bound of alpha for input_real_shu


  // For NAM data, use the following command and initial values...
  // does NOT work very well as of 04/30/2023
  //a.out input1_nam_spearman outpara1_nam_spearman  1.0e-3   1.0
  a1=0.0000001 ;  // lower bound of lambda for input_nam
  b1=10.0 ;     // upper bound of lambda for input_nam

  a2=0.0000001 ;  // lower bound of alpha for input_nam
  b2=10.0 ;     // upper bound of alpha for input_nam


  //--------------------------------------------------------------------
  //--------------------------------------------------------------------

  // count the number of data
  n=0 ; 
  while(fscanf(fi, "%d %lf", &val, &valy)!=EOF)  
    n++ ;
  rewind(fi) ;

  // read input data
  dx=(int *)malloc(sizeof(int)*n) ;
  dy=(double *)malloc(sizeof(double)*n) ;
  for(i=0 ; i < n ; i++) {
    fscanf(fi, "%d %lf", dx+i,  dy+i) ;
    //    printf("%d %lf\n", dx[i], dy[i]);
  }
  fclose(fi) ;

  // alpha=0.1 ;   // initial arbitrary value of alpha
  // alpha=0.3 ;   // initial arbitrary value of alpha
  count=0 ;   // for convergence
  olambda=oalpha=-1.0 ;   // inital values
  //  for(i=0 ; i < 300 ; i ++) {
  for(i=0 ; i < 800 ; i ++) {
    // find root lamdba for given alpha
    lambda=rtbis(pfun1, a1, b1, alpha, acc, dx, dy, n) ;
      
    alpha=rtbis(pfun2, a2, b2, lambda, acc, dx, dy, n) ;

    printf("%lf  %lf\n", lambda, alpha) ;

    if(olambda==lambda && oalpha==alpha)
      count++ ;
      
    if(count==5) {
      fprintf(fo, "%lf %lf\t", lambda, alpha) ;
      break ;
    }

    olambda=lambda ;
    oalpha=alpha ;
  }

  //if(i >= 300) { printf("iteration is not enough...\n") ;  exit(1) ; }
  if(i >= 800) { printf("iteration is not enough...\n") ;  exit(1) ; }

  fprintf(fo, "%lf %lf\n", lambda, alpha) ;

  free(dx) ;
  free(dy) ;
  return 0 ;
}


/*                                                                             
 *                                                                             
 * Using bisection, find the root of a function "func" known to lie between x1
 * and x2. The root, returned as "rtbis", will be refined until it accuracy is
 * \pm x acc 
 * 
 */
double rtbis(double (*func)(int*, double*, int, double, double), double x1, double x2, double para, double xacc, int *dax, double *day, int n)
{
  int j ;
  double dx, f, fmid, xmid, rtb ;

  f=(*func)(dax, day, n, x1, para) ;
  fmid=(*func)(dax, day, n, x2, para) ;

  if(f*fmid >= 0.0) {
    printf("root must be bracketed for bisection in rtbis %lf\n", para) ;
    printf("%lf %lf\n", f, fmid) ;
    exit(1) ;
  }

  rtb = f < 0.0 ? (dx=x2-x1, x1) : (dx=x1-x2, x2) ;

  for(j=1 ; j <= JMAX ; j++) {
    //    printf("%2d %lf %lf", j, rtb, rtb+dx) ;

    fmid=(*func)(dax, day, n, xmid=rtb+(dx *= 0.5), para) ;
    if(fmid <= 0.0)
      rtb=xmid ;

    //    printf(" %lf %10lf\n", xmid, fmid) ;

    if(fabs(dx) < xacc || fmid == 0.0 )
      return rtb ;
  }
  printf("Too many bisection in rtbis\n") ;
  exit(1) ;
  
  return 0.0 ;
}
