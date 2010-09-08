#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #define MATHLIB_STANDALONE 1 */
/* #include <Rmath.h> */
#include "mat_vec.h"
#include "array.h"
#include "order.h"
#include "quantile.h"
#define Inf 1e+140
#define SQ(x) ((x)*(x))


void spherize(double **x, int nrow, int ncol);

int kmeandirstarts(int n, int p, int nclus, double **x, double **Mu, 
		   int reqsvd, int startmeth, int nrandoms);

void meandispersion(double **x, int n, int p, double *mu, double *ltsigma);

int initials(double **x,int n,int p,int nclass,int *nc,
		 double **Mu,double **LTSigma,int *class);

void kmean_dirs(double **a, int m, int n, double **c, int k, int *ic1,int *nc,
	   int iter,double *wss,int *ifault,double *normc);

int run_kmndirs(double **X, int n, int p, double **mu, int K, int *id, 
		int iter,double *normc, int *nc, int nrand)
{
  int i=0;

  if(K == 1){
    double **ltsigma;
    for(i = 0; i < n; i++) id[i] = 0;
    
    MAKE_MATRIX(ltsigma, K, p*(p+1)/2);
    
    i = initials(X, n, p, K, nc, mu, ltsigma, id);
    spherize(mu, K, p);

    FREE_MATRIX(ltsigma);
  }
  else{
    int IFAULT, reqsvd = 0, startmeth = 0;
    double *WSS;
    
    /* 1. Get starting values.  */
    i = kmeandirstarts(n, p, K, X, mu, reqsvd,startmeth,nrand);

    /*2. Run k-mean directions  */
    MAKE_VECTOR(WSS,K);
/*    MAKE_VECTOR(nc,K);*/
    
    kmean_dirs(X, n, p, mu, K, id, nc,iter, WSS,(&IFAULT),normc);

/*    FREE_VECTOR(nc);*/
    FREE_VECTOR(WSS);
}

  return i;
}

