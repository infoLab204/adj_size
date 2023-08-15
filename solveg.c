/*                                                                            
 * find a solution to a given equation
 *
 * $a.out  para3_input_spearman  1.0e-5 > gout3_spearman  // model-free data
 * $a.out  para1_input_spearman  1.0e-5 > gout1_spearman  // real dat
 *
 * (note) para3_input_spearman consists of dependent alpha, lambda, and independent alpha, lambda
 *  
 */

#define JMAX 40   // maximum allowed number of bisections  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double ffunc(double alpha, double lambda, double x, double con) ;
double rtbis(double (*func)(double, double, double, double), double a2, double b2, double x1, double x2, double con, double xacc) ;

int main(int argc, char **argv)
{
  FILE *fi ;
  double a1, a2, b1, b2, x1, x2, x, acc, val1, root ;
  double (*pfun)(double a, double b, double x, double val) ;

  if(argc!=3) {
    printf("enter correct parameters...\n") ;
    exit(1) ;
  }

  // a1, b1: alpha and lambda of correlated (or real) case
  // a2, b2: alpha and lambda of uncorrelated (random) case 
  // (eg) a1, b1: alpha and lambda of para3_real_pearson
  // (eg) a2, b2: alpha and lambda of para3_random_pearson
  fi=fopen(argv[1], "r") ;   // read input parameters
  fscanf(fi, "%lf %lf %lf %lf", &a1, &b1, &a2, &b2);
  acc=atof(argv[2])  ;  // error or accuracy of the bisection method   

  pfun=ffunc ; // function pointer

  for(x=10 ; x <= 850 ; x+=2) {   // sample size of correlated or real model
    // evaluate function of correlated or real model
    // (note) the value of 0.0 
    val1=ffunc(a1, b1, x, 0.0); // 'x' is the numeber of samples...
 
    x1=1.0 ;      // lower bound of sample size
    x2=1000.0 ;   // upper bound of sample size

    // find root of the corresponding sample size of uncorrelated model
    root=rtbis(pfun, a2, b2, x1, x2, val1, acc) ;
    printf("%lf\t%lf\n", x, root) ;  // root is the corresponding sample size
  }
  return 0 ;
}

////////////////////////////////
// evaluate the model function
double ffunc(double alpha, double lambda, double x, double con)
{
  double pp, ff ;

  pp=pow(x, 1.0-alpha) ;
  ff=lambda/(1.0-alpha)*(1.0-pp)+log(lambda)-alpha*log(1.0*x)-con ;

  return ff ;
}
///////////////////////////////

/*                                                                             
 *                                                                             
 * Using bisection, find the root of a function "func" known to lie between x1
 * and x2. The root, returned as "rtbis", will be refined until it accuracy is
 * \pm x acc 
 * 
 */
double rtbis(double (*func)(double, double, double, double), double a2, double b2, double x1, double x2, double con, double xacc)
{
  int j ;
  double dx, f, fmid, xmid, rtb ;

  f=(*func)(a2, b2, x1, con) ;      // x1 is lower limit
  fmid=(*func)(a2, b2, x2, con) ;   // x2 is upper limit

  if(f*fmid >= 0.0) {
    printf("root must be bracketed for bisection in rtbis...\n") ;
    printf("%lf %lf\n", f, fmid) ;
    exit(1) ;
  }

  rtb = f < 0.0 ? (dx=x2-x1, x1) : (dx=x1-x2, x2) ;

  for(j=1 ; j <= JMAX ; j++) {
    //    printf("%2d %lf %lf", j, rtb, rtb+dx) ;
    fmid=(*func)(a2, b2, xmid=rtb+(dx *= 0.5), con) ;
    if(fmid <= 0.0) rtb=xmid ;
    //    printf(" %lf %10lf\n", xmid, fmid) ;
    if(fabs(dx) < xacc || fmid==0.0) 
      return rtb ;
  }
  printf("Too many bisection in rtbis\n") ;
  exit(1) ;
}
