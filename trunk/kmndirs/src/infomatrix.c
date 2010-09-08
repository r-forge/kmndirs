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


int svdd(double **a, int m, int n, double *d, double **u, double **v);
int print_nc_min(int *id, int n, int K);
int sphere_proj(double **X, double **xprime, int n, int m, double **mu, double **muprime, int K);
double median(int n, double *x);

double obj_function(double **X, double **mu, int *id, int n, int p)
{
  int i,j;
  double sum=0.,temp;

  for(i = 0; i < n; i++){
    temp=0.;
    for(j = 0; j < p; j++)  temp += SQ(X[i][j]-mu[id[i]][j]);
    /*    sum += sqrt(temp);*/
    sum += temp;
  }
  return sum;
}

int  sphere_proj(double **X, double **xprime, int n, int m, double **mu, double **muprime, int K)
{
      double **U, **V, *D, **xsamp;
      int i,j, d2, nm,ll;

      if(m < n) nm = m;
      else nm = n;

      MAKE_MATRIX(xsamp, nm, m);
      MAKE_MATRIX(U, nm, m);
      MAKE_MATRIX(V, m, m);
      MAKE_VECTOR(D, m);

      for(i = 0; i < nm; i++){
	for(j = 0; j < m; j++)	xsamp[i][j]=X[i][j];
      }

      d2 = svdd(xsamp, nm, m, D, U, V);

      FREE_MATRIX(xsamp);
      FREE_MATRIX(U);
      FREE_VECTOR(D);


/*      for(ll = 0; ll < m; ll++) printf("%lf ",(X[0][ll] * V[ll][0]));*/
      
      for(i = 0; i < n; i++){
	for(j = 0; j < (m - 1); j++){
	  xprime[i][j]=0.;
	  for(ll = 0; ll < m; ll++) xprime[i][j] += (X[i][ll] * V[ll][j]);
	}
      }
      for(i = 0; i < K; i++){
	for(j = 0; j < (m - 1); j++){
	  muprime[i][j]=0.;
	  for(ll = 0; ll < m; ll++) muprime[i][j] += mu[i][ll] * V[ll][j];
	}
      }

      FREE_MATRIX(V);

      return i;
}



int print_nc_min(int *id, int n, int K){
  int i, *nc, min_nc=n+1;

  MAKE_VECTOR(nc,K);

  for(i=0;i<K;i++) nc[i]=0;

  for(i=0;i<n;i++) nc[id[i]]++;

  for(i=0;i<K;i++){
    if(min_nc > nc[i]) min_nc=nc[i];
  }

  FREE_VECTOR(nc);

  return min_nc;
}



