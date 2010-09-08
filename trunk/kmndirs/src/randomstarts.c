#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/* #define MATHLIB_STANDALONE 1 */
/* #include <Rmath.h> */
#include "array.h"
#include "mat_vec.h"
#define Inf 1e+140
#define SQ(x) ((x)*(x))

double dotproduct(double *x, double *y, int n);
double obj_function(double **X, double **mu, int *id, int n, int p);


int assign_closest(double *X,int p,int nclass,double **Mu);


int srswor(int n, int k, int *y)
{
    /* Provide k out of n indices sampled at random without replacement 
       Author: xxxxxxxxxx
    */
  
    if(k > n) {
	REprintf("Error: k greater than n in srswor()");
	return 1;
    }
    else {
	
	int i, j;
	int *x;
	
	MAKE_VECTOR(x, n);
	for (i = 0; i < n; i++)	x[i] = i;

	GetRNGstate();
    
	for (i = 0; i < k; i++) {
	    j = n * unif_rand();
	    /*      printf("%d ",j);*/
	    y[i] = x[j];
	    x[j] = x[--n];
	}

	PutRNGstate();

	/*    printf("\n\n\n");*/

	FREE_VECTOR(x);
    }
    return 0;
}


int rand_inits(double **X, int n, int p, int K, double **MU, int *id, 
	       int nstarts, double *lval)
{
  /*
    X: data
    n,p: obs,dims
    K: # clusters
    MU: Kxp matrix of cluster means
    id: classifications
    E: minimum cluster representation
    SBI: Indicator whether SBI should be used instead of mean resultant
    kmeth: method of estimating kappa
    maxtries: # of tries to obtain keep==K
  */

  int flag=0,i,j,k,ns,*samp,dumm,*nc, *optid;
  double like, maxlike, **optMu;


  MAKE_VECTOR(samp,K);
  MAKE_MATRIX(optMu,K,p);
  MAKE_VECTOR(optid, n);
  MAKE_VECTOR(nc,K);



  maxlike = Inf;

  /*  for(ns=0;(ns<nstarts)&&(flag==0);ns++){*/
  ns=0;
  while(ns<nstarts){
  /* First get random sample with required min representation.  */
    
    dumm=srswor(n,K,samp);


    for(i=0;i<K;i++){
      for(j=0;j<p;j++)  MU[i][j]=X[samp[i]][j];
    }


    /* Next step is to make one run of the algorithm to determine "optimal"
       starting value.*/
    for(k=0;k<K;k++) nc[k]=0;
    for(i=0;i<n;i++){
      id[i] = assign_closest(X[i], p, K, MU);
      nc[id[i]]++;
    }
    
    for(k=0;(k<K)&&(flag==0);k++){
      if(nc[k]==0){
        flag=1;
/*        printf("Flagged! Not good! R = %d\n",ns);*/
      }
    }

    if(flag==0){
      ns++;
      like = 0.;
      
      like =  obj_function(X, MU, id, n,p);
      
      if(like<maxlike){
	for(k=0;k<K;k++){
	  for(j=0;j<p;j++){
	    optMu[k][j]=MU[k][j];
	  }
	}
	for(i=0;i<n;i++) optid[i]=id[i];
	maxlike=like;
      }
    }

    if(flag==1){
/*      printf("Obtaining new random.seed\n");*/
      FILE *fran;
      unsigned int seed1, seed2;
      fran = fopen("random.seed","r");
      fscanf(fran, "%u %u", &seed1, &seed2);
      seed1++;
      seed2++;
      fclose(fran);
      
      fran=fopen("random.seed", "w");
      fprintf(fran, "%d %d\n", seed1, seed2);
      fclose(fran);
      flag=0;
    }
  }

 /* if(flag==0){*/
    for(k=0;k<K;k++){
      for(j=0;j<p;j++){
	MU[k][j]=optMu[k][j];
      }
    }
    for(i=0;i<n;i++) id[i]=optid[i];
    (*lval) = maxlike;
/*  }*/

  FREE_MATRIX(optMu);
  FREE_VECTOR(samp);
  FREE_VECTOR(nc);
  FREE_VECTOR(optid);

  REprintf("Finished Random Starts\n");

  return flag;
}

