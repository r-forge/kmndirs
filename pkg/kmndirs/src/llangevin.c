#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "array.h"
/* #define MATHLIB_STANDALONE 1 */
/*It is essential to have this before the call 
  to the Rmath's header file because this decides
  the definitions to be set. */
#include <Rmath.h>
#define PI 3.141593
#define Inf 1e+140
#define SQ(x) ((x) * (x))

void standardize(double **x, int nrow, int ncol);


void make_unit(double **x, int n, int p)
{
/*
  This function ensures the n rows of x lie on the unit p-sphere
*/
	int i, j;
	double sums2;

	for(i = 0; i < n; i++) {
		sums2 = 0.;
		for(j = 0; j < p; j++) sums2 += SQ(x[i][j]);
		for(j = 0; j < p; j++) x[i][j] /= sqrt(sums2);
	}
	
	return;
}


double besselI(double x, double ord, double off){
/* 
   Calls the Rmath library's bessel_i  function for a given value and order.
   x: value
   ord: order
   off: set to 2 if returning exp(-x)I_ord(x) ...recommended.
   
   Note that although some approximation are placed in this function, it
   is not very stable for large dimensions.
*/

  double bi;

  if(x < 1e-5) x = 1e-4;


  /*  if(off>1.){
    if(x>1000.)     bi=1./sqrt((2.*x*PI));
    else      bi = bessel_i(x, ord, off);
  }
  else{
    if(x>1000.)      bi = exp(x) / sqrt(2.*x*PI);
    else      bi = bessel_i(x, ord, off);
  }
  */

  bi = bessel_i(x, ord, off);

  return bi;
}

double Iaprx(double k,double p,int logged){
  double bess=0.,value;

  value = 1. - (4.*p*p-1.)/(8.*k) 
    + (4.*p*p-1.)*(4.*p*p-9.)/(128.*k*k) - (4.*p*p-1)*(4.*p*p-9.)
    *(4.*p*p-25.)/(1024.*k*k*k);

  if( value > 0 ){
    if(!logged)    bess = (1./sqrt(2.*PI*k)) * exp(k) * value;
    if(logged)     bess = -.5*log(2.*PI*k)+k+log(value);
  }
  else{
    if(k <= p){
      if(!logged)  bess = besselI(k, p,1);
      if(logged)   bess = log(besselI(k, p,1));
    }
    else{
      if(!logged)  bess = besselI(k, p,2)*exp(k);
      if(logged)   bess = log(besselI(k, p,2)) + k;
    }
  }

/*  printf("besselI = %lf\n",bess);*/

  return bess;
}





double adk(double p, double kappa){
/*
  Compute the A_d(k) value which is a ratio of bessel functions.
*/
  double num,dem,A;

  num = Iaprx(kappa,p/2.,1);
  dem = Iaprx(kappa,p/2-1,1);
  A = exp(num - dem);

  /*
  if(kappa<1500){
    num = besselI(kappa,p/2,2.);
    dem = besselI(kappa,p/2-1,2.);
    
  }
  else{
    num=1.;
    dem=1.;
  }

  A=num/dem;
*/
  return A;
}

double dunit_sphere(int p)
{
  double ptemp,gam,sphere;
  gam = 1.0;
  ptemp = p / 2. + 1.;
  while(ptemp > 1.0){
    ptemp -= 1.0;
    gam *= ptemp;
  }
  if(ptemp != 1.0) gam *= 0.5 * sqrt(PI);

  sphere = -0.5 * p * log(2. * PI) + log(gam);

  return sphere;
}

double logdlang(double *x, double *mu, double kappa, double constant, int p)
{
/*
  Internal function used to calculate the "kernal" of the Langevin.
*/
  int j;
  double dl;

  dl = constant;

  for(j = 0; j < p; j++){
    dl += kappa * x[j] * mu[j];
  }

  return dl;
}


double llikecoL(double *x, double *nu, double kappa, int p)
{

  int j,d = p - 1;
  double dl;

  dl = -1.*(d/2.) * log(2.*PI) - (1.- d/2.)*log(kappa) - 
    Iaprx(kappa,(d/2.-1.),1);


  for(j = 0; j < p; j++){
    dl += kappa * x[j] * nu[j];
  }

  return dl;
} 
