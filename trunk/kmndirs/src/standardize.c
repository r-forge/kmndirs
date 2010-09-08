#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define BIG 1e+40
#include "array.h"

/* This function is used to force a vector to be a unit vector with mean zero.
i.e sum(m)=0 and sum(m^2)=1 
  IMPORTANT! THIS CHANGES THE VALUES IN THE INPUT VECTOR. 
*/


void unitize(double *m,int lengthm)
{

  int i;
  double sum;
  double mean=0.0;
  sum=0.0;
  for(i=0;i<lengthm;i++){
    mean+=m[i];
  }

  mean/=1.0*lengthm;

  for(i=0;i<lengthm;i++){
    m[i]-=mean;
    sum+=m[i]*m[i];
  }

  for(i=0;i<lengthm;i++){
    m[i]/=sqrt(sum);
  }
  
  return;
}

void standardize(double **x, int nrow, int ncol)
{

  /* This function standardizes the rows of x so that they have zero mean 
     and lie on the unit ncol-dimensional sphere */

  double xmu, xvar;
  int i, j;

  for (i=0; i<nrow; i++) {
    xmu = 0.0;
    for (j=0; j<ncol; j++) xmu += x[i][j]/ncol;
    for (j=0; j<ncol; j++) x[i][j] -= xmu;
    xvar = 0.0;
    for (j=0; j<ncol; j++) xvar += x[i][j] * x[i][j];
    if (xvar > 0) {
      for (j=0; j<ncol; j++) x[i][j] /= sqrt(xvar);
    }
  }

  return;
}


