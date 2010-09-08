#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/* #define MATHLIB_STANDALONE 1 */
/* #include <Rmath.h> */
#include <string.h>
#include "array.h"
#include "mat_vec.h"

/* #define PI 3.141593 */
#define Inf 1e+143
#define SQ(x) ((x) * (x))


int quantile(int n,double *x,double *p,double *q, int numqs);

double dotproduct(double *x, double *y, int n)
{
  int i;
  double sum=0.;
  for(i = 0; i < n; i++) sum += x[i]*y[i];
  return sum;
}

double median(int n, double *x){
/* Calls the quantile function to get a sample median. */
  int dum;
  double med,p=0.50;
  dum=quantile(n,x,&p,&med,1);
  return med;
}

