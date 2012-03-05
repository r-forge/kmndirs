#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "array.h"

/* #define PI 3.141593 */
#define Inf 1e+140

double chisqstatlt(int p, double *X,double *mu,double *ltSigma,double *detSig);

int eigens(double *A, double *EV, double *E, int n);

double dlmvnorm(double *x, int p, double *mu, double *LTsigma)
{
  /* Calculates the p-variate Gaussian log-density at the p-variate point x
     with mean mu and the lower triangle (px(p+1)/2) LTsigma of the 
     pxp-dimensional dispersion matrix Sigma. */
  double detsig,exponent;
  exponent=-chisqstatlt(p,x,mu,LTsigma,&detsig);
  if (detsig>0) {
    exponent*=0.5;
    exponent-=0.5*log(detsig)+0.5*p*log(2*PI);
  }
  return(exponent);
}

double logdeterminant(double *temp,int n)
/* This function calculates the log-determinant of symmetric (hence 
   nonnegative-definite) matrices in packed storage*/
{
  int i;
  double *eigvec,*eigval,sump=0.;
  MAKE_VECTOR(eigvec,n*n);
  MAKE_VECTOR(eigval,n);
  if (!eigens(temp,eigvec,eigval,n)) {
    for(i=0;(i<n) && (eigval[i]>0);i++) sump+=log(eigval[i]);
    if (i==0) sump=-Inf;
  }
  else REprintf("Serious error in log-determinant calculations\n");
  FREE_VECTOR(eigval);
  FREE_VECTOR(eigvec);
  return sump;
}  

double dlmvnorm_singular(double *x, int p, double *mu, double *LTsigma)
{
  /* Calculates the p-variate Gaussian log-density at the p-variate point x
     with mean mu and the lower triangle (px(p+1)/2) LTsigma of the 
     pxp-dimensional dispersion matrix Sigma. 
     Note that this function actually incorporates the case for singular 
     determinants
  */
  double *eivec,*eival,value=0;
  int i,ind=0;

  MAKE_VECTOR(eivec,p*p);
  MAKE_VECTOR(eival,p);
  i=eigens(LTsigma,eivec,eival,p);  
  if (eival[0]==0) {      /* only possible if LTsigma is all zero for
			     symmetric matrices, which means 
			     that the distribution is degenerate at mu */
    for(i=0;((i<p) && (!ind));i++) if (x[i]!=mu[i]) ind =1;
    if (ind) value=-Inf;
    else value=0;
  }
  else {
    int j,dmin;
    double *y,*z,sum=0,sump=0;
    for(i=0;((i<p) &&  (eival[i]>0));i++) sum+=eival[i];
    for(i=0;((i<p) && (sump<1.) &&(eival[i]>0));i++) {
      sump+=eival[i]/sum;
      value-=0.5*log(eival[i]);
    }
    dmin=i;
    MAKE_VECTOR(y,p);
    for(i=0;i<p;i++) y[i]=x[i]-mu[i];
    MAKE_VECTOR(z,dmin);
    for(i=0;i<dmin;i++) z[i]=0;
    for(i=0;i<dmin;i++)  {
      for(j=0;j<p;j++) z[i]+=eivec[j*p+i]*y[j];
    }
    FREE_VECTOR(y);
    for(i=0;i<dmin;i++) value-=0.5*z[i]*z[i]/eival[i];
    FREE_VECTOR(z);
    value-=0.5*dmin*log(2*PI);
  }
  FREE_VECTOR(eival);
  FREE_VECTOR(eivec);
  return value;
}

double mixllhd(int p,int k,double *x,double *pi,double **Mu,double **LTSigma) 
{
  /* This function calculates the likelihood of a mixture of p-variate 
     Gaussians at a p-variate point X. Each density has prior class 
     probabilities given by pi, and has mean given by the kxp-matrix Mu 
     and the kxp(p+1)/2 matrix of corresponding lower triangular dispersion 
     matrix vectors */
  int i;
  double temp,sum=0.;
  for (i=0;i<k;i++) {
    temp=dlmvnorm(x,p,Mu[i],LTSigma[i]);
    /* temp=dlmvnorm_singular(x,p,Mu[i],LTSigma[i]);*/
    sum+=pi[i]*exp(temp);
  }
  return sum;
}

double lnlikelihood(int n,int p,int k,double *pi,double **X,double **Mu,
                    double **LTSigma)
{
  int i;
  double llhd,temp;
  llhd=0.;
  for(i=0;i<n;i++) {
    temp=mixllhd(p,k,X[i],pi,Mu,LTSigma);
    llhd+=log(temp);
  }
  return(llhd);
}


