/* using an extension of Maitra's~(2009) algorithm for sphered data 

  modified by xxxx, xxxx and xxxx
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mat_vec.h"
#include "array.h"
#include "order.h"
#include "quantile.h"
#define Inf 1e+140
#define SQ(x) ((x)*(x))
#define PI 3.141593

void hclassify(int n,int m, double **x,int hcrit,int nclass,int *class);

void kmeans(double **a, int m, int n, double **c, int k, int *ic1,int *nc,
	    int iter,double *wss,int *ifault);

int sort(int n,double *x);
int svdd(double **a, int m, int n, double *d, double **u, double **v);
void unitize(double *m,int lengthm);
int polaroid(int p, double *x, double *R, double *theta);
int assign_closest(double *X, int p, int nclass, double **Mu);
double pposymatdet(int N,double *A, char UPLO);
void spherize(double **x, int nrow, int ncol);

int rand_inits(double **X, int n, int p, int K, double **MU, int *id,
	       int nstarts, double *lval);
double dotproduct(double *x, double *y, int n);
double obj_function(double **X, double **mu, int *id, int n, int p);


double determinant(double *LTSigma,int n)
{
/*  double dum;
  pposymatinv(n,LTSigma,'U',&dum);
  return dum;*/

  return pposymatdet(n,LTSigma,'U');

}

void meandispersion(double **x, int n, int p, double *mu, double *ltsigma)
{
  /* This routine calculates the mean and dispersion of a homogeneous sample */
  int i,j,l;

  for (i=0;i<(p*(p+1)/2);i++) ltsigma[i]=0.;
  for (i=0;i<p;i++) mu[i]=0.;
  for (i=0;i<n;i++) {
    for (j=0;j<p;j++) mu[j]+=x[i][j];
  }
  for (j=0;j<p;j++) mu[j]/=n;
  for (i=0;i<n;i++) {
    for (j=0;j<p;j++) {
      for (l=0;l<=j;l++) ltsigma[j*(j+1)/2+l]+=(x[i][j]-mu[j])*(x[i][l]-mu[l]);
    }
  }
  if(n>1) {
    for (j=0;j<(p*(p+1)/2);j++) ltsigma[j]/=n-1;
  }
  return;
}

void meandispersion_dir(double **x, int n, int p, double *mu, double *ltsigma)
{
  /* This routine calculates the mean and dispersion of a homogeneous sample */
  int i,j,l;

  for (i=0;i<(p*(p+1)/2);i++) ltsigma[i]=0.;
  for (i=0;i<p;i++) mu[i]=0.;
  for (i=0;i<n;i++) {
    for (j=0;j<p;j++) mu[j]+=x[i][j];
  }
  for (j=0;j<p;j++) mu[j]/=n;

  unitize(mu, p);  /* Enforce constraint!  */

  for (i=0;i<n;i++) {
    for (j=0;j<p;j++) {
      for (l=0;l<=j;l++) ltsigma[j*(j+1)/2+l]+=(x[i][j]-mu[j])*(x[i][l]-mu[l]);
    }
  }
  if(n>1) {
    for (j=0;j<(p*(p+1)/2);j++) ltsigma[j]/=n-1;
  }
  return;
}


int initials(double **x,int n,int p,int nclass,int *nc,
	      double **Mu,double **LTSigma,int *class)
{
  double **y;
  int i,j,k,l,m=1;
  for (i=0;i<nclass;i++) {
       nc[i]=0;
       for(l=0;l<n;l++) if (class[l]==i) nc[i]++;
  }
  for(i=0;i<nclass;i++) {
    if(nc[i]>p) m*=1;
    else m*=0;

    MAKE_MATRIX(y,nc[i],p);
    k=0;
    for(l=0;l<n;l++) {
      if (class[l]==i) {
	for (j=0;j<p;j++)   y[k][j]=x[l][j];
	k++;
      }
    }
    meandispersion(y,nc[i],p,Mu[i],LTSigma[i]);
    FREE_MATRIX(y);
  }
  return m;
}


int initials_dir(double **x,int n,int p,int nclass,int *nc,
	      double **Mu,double **LTSigma,int *class)
{
  double **y;
  int i,j,k,l,m=1;
  for (i=0;i<nclass;i++) {
       nc[i]=0;
       for(l=0;l<n;l++) if (class[l]==i) nc[i]++;
  }
  for(i=0;i<nclass;i++) {
    if(nc[i]>p) m*=1;
    else m*=0;

    MAKE_MATRIX(y,nc[i],p);
    k=0;
    for(l=0;l<n;l++) {
      if (class[l]==i) {
	for (j=0;j<p;j++)   y[k][j]=x[l][j];
	k++;
      }
    }
    meandispersion_dir(y,nc[i],p,Mu[i],LTSigma[i]);
    FREE_MATRIX(y);
  }
  return m;
}

/* write the digits of n in its broken down "base" into buf */
void break_down(int n, int *base, int *buf, int buflen)
{
        int i;
        for (i = 0; i < buflen; i++) {
	  buf[i] = n % base[i];
	  n /= base[i];
        }
}


int assign_closest(double *X, int p, int nclass, double **Mu)
{
	int j, l, class = 0;
	double temp, dum1;
	
	temp = Inf;
	for (l = 0; l < nclass; l++) {
		dum1 = 0.;
		for(j = 0; j < p; j++) dum1 += SQ(X[j] - Mu[l][j]);
		if (dum1 < temp) {
			temp = dum1;
			class = l;
		}
	}
	return class;
}

double **eliminulls(double **x, int n, int p, int *nclass, double **Mu, int kk)
{ 
  /*This routine eliminates all those centers which have representation less
    than or equal to kk and returns the remaining centers */

	int i, j, *nc, k=0, newtotcl = (*nclass);
	double **Centers;

	MAKE_VECTOR(nc, (*nclass));
	for(i = 0; i < (*nclass); i++) nc[i] = 0;
	for(i = 0; i < n; i++) nc[ assign_closest(x[i], p, (*nclass), Mu)]++;
	for(i = 0; i < (*nclass); i++) if(nc[i] <= kk) newtotcl--;
	
	MAKE_MATRIX(Centers, newtotcl, p);
	for(i = 0; i < (*nclass); i++){
		if (nc[i] > kk) {
			for(j = 0; j < p; j++) Centers[k][j] = Mu[i][j];
			k++;
		}
	}
	(*nclass) = newtotcl;
	FREE_VECTOR(nc);
	return Centers;
}

void minmax(int n, double *x, double *minmax) /* finds the mininmum and 
						 maximum of x */
{
	int i;
	
	minmax[0] = x[0];
	minmax[1] = x[0];
	for (i = 1; i < n; i++) {
		if (x[i] < minmax[0]) minmax[0] = x[i];
		else {
			if (x[i] > minmax[1]) minmax[1] = x[i];
		}
	}
	
	return;
}

void partialfill(double **full, double **partial, int nf, int np, int m){
  int i, j;
  for(i = 0; i < np; i++){
    for(j = 0; j < m; j++){
      full[i][j] = partial[i][j];
    }
  }
  return;
}


int get_polars_range_for_last_angl(int n, int p, double **x, double **theta,
				   double *psi)
{
	int i, j, whichangl = 0; 
	/* whichangl = 0 if the data in enclosed by the 
	   outer limits, 1 if the data is outside the 
	   limits. Note that only the last angle has 
	   this issue, being on the circle: the others
	   are on portions of half-circles. */			   
	
	double r, *tmp, jump;

	MAKE_VECTOR(tmp, n);

	for (i = 0; i < n; i++) j = polaroid(p, x[i], &r, theta[i]);
	for (i = 0; i < n; i++) tmp[i] = theta[i][p-2];

	sort(n, tmp);
	
	psi[0] = tmp[0];
	psi[1] = tmp[n-1];
	jump = 2 * PI + tmp[0] - tmp[n-1];
	whichangl = 0;
	
	for (i = 0; i < (n - 1); i++) { 
		if ( (tmp[i+1] - tmp[i]) > jump) {
			jump = tmp[i+1] - tmp[i];
			psi[0] = tmp[i];
			psi[1] = tmp[i+1];
			whichangl = 1;
		}
	}
	
	FREE_VECTOR(tmp);

	return whichangl;
}


int starts_in_polar_domain(int n,int m,double **x,int nclus,int *ningrp,
			   int *grpids)
{    
  double **cent, *dum1, *probs, *qtl, *cctr, **dumy, **dumcctr,
    *wss1, **ccent, **mu, **ltsigma, **y, **tempcent;
  int i, j, k, *counts, *ncl, totcl = 1, sumcl = 0, *dum2, *buf, *ning,
    kopt, *grclass, *cum1, maxiter = 10, ind = 1,complete = 0, 
    index = 0,box,partn;   
  
  MAKE_VECTOR(ncl, m);
  i = (int) ceil( pow((m*nclus), 1.0 / m) );
  for(j = 0; j < m; j++) ncl[j] = i;
  for(j = 0; j < m; j++) {
    totcl *= ncl[j];
    sumcl += ncl[j];
  }
	
  MAKE_MATRIX(y, n, m);
  MAKE_VECTOR(cctr, m);
  MAKE_VECTOR(dum1, m*(m+1)/2);
  meandispersion(x, n, m, cctr, dum1);
  
  for(i = 0; i < n; i++) {
    for(j = 0; j < m; j++) y[i][j] = (x[i][j] - cctr[j]) / sqrt(dum1[j*(j+1)/2+j]);
  }
  FREE_VECTOR(dum1);
  FREE_VECTOR(cctr);
  
  MAKE_VECTOR(dum1,n);
  MAKE_VECTOR(cctr,sumcl);
  MAKE_MATRIX(dumy,n,1);
  MAKE_VECTOR(dum2,n);
  
  k = 0;
  for(i = 0; i < m; i++) {
    MAKE_VECTOR(ning, ncl[i]);
    MAKE_VECTOR(qtl, ncl[i]);
    MAKE_MATRIX(dumcctr, ncl[i], 1);
    MAKE_VECTOR(wss1, ncl[i]);
    for(j = 0; j < n; j++) dum1[j] = y[j][i];
    MAKE_VECTOR(probs, ncl[i]);
    for(j = 0; j < ncl[i]; j++) probs[j] = (2*j+1)/(2.*ncl[i]);
    j = quantile(n, dum1, probs, qtl, ncl[i]);
    for(j = 0; j < n; j++) dumy[j][0] = y[j][i];
    for(j = 0; j < ncl[i]; j++) dumcctr[j][0] = qtl[j];
    kmeans(dumy, n, 1, dumcctr, ncl[i], dum2, ning, maxiter, wss1,
	   &kopt);
    for(j=0;j<ncl[i];j++) cctr[k++]=dumcctr[j][0];
    FREE_VECTOR(probs);
    FREE_VECTOR(wss1);
    FREE_VECTOR(ning);
    FREE_VECTOR(qtl);
    FREE_MATRIX(dumcctr);
  }
  FREE_VECTOR(dum2);
  FREE_MATRIX(dumy);
  FREE_VECTOR(dum1);
	
  MAKE_VECTOR(cum1,m);
  cum1[0]=0;
  for(i=1;i<m;i++) cum1[i]=cum1[i-1]+ncl[i-1];



  box=n*((int)ceil(sqrt(m)));

  MAKE_MATRIX(ccent,box,m);
	
  MAKE_VECTOR(buf,m);

  while(complete <= totcl){
    for(i = index; i < box; i++){
      complete++;
      break_down(i, ncl, buf, m);
      for(j=0;j<m;j++) ccent[i][j]=cctr[cum1[j]+buf[j]];
    }
    partn = box;
    tempcent = eliminulls(y, n, m, &partn, ccent,0);
    partialfill(ccent, tempcent, box, partn, m);
    FREE_MATRIX(tempcent);
    index = partn;
  }

  FREE_VECTOR(cum1);
  FREE_VECTOR(cctr);
  FREE_VECTOR(buf);
  FREE_VECTOR(ncl);
  
  cent=eliminulls(y,n,m,&partn,ccent,0);
  totcl = partn;

  FREE_MATRIX(ccent);
  
  if(totcl>=nclus) {
    
    ind=0;
    MAKE_VECTOR(counts,totcl);
    MAKE_VECTOR(wss1,totcl);
    kmeans(y,n,m,cent,totcl,grpids,counts,maxiter,wss1,&k);
    FREE_VECTOR(wss1);
    FREE_VECTOR(counts);
    
    MAKE_VECTOR(grclass,totcl);
    hclassify(totcl,m,cent,2,nclus,grclass);
    
    MAKE_MATRIX(mu,totcl,m);
    MAKE_MATRIX(ltsigma,totcl,m*(m+1)/2);
    i=initials(cent,totcl,m,nclus,ningrp,mu,ltsigma,grclass);
    FREE_VECTOR(grclass);
    FREE_MATRIX(cent);
    for(i=0;i<n;i++) grpids[i]=assign_closest(y[i],m,nclus,mu);
    FREE_MATRIX(mu);
    FREE_MATRIX(ltsigma);
  }
  else FREE_MATRIX(cent);
  FREE_MATRIX(y);
  return ind;
}




int jumpstarters(int n, int m, double **y, int nclus, int *ningrp, int *grpids,
		 int reqsvd, int startmeth, int nrandoms)
{

  int i,ind = 0, *id;
  double maxval = Inf, temp, **mu, **ltsigma;

  MAKE_MATRIX(mu,nclus,m);
  MAKE_VECTOR(id,n);

  if(nrandoms > 0){
    ind = rand_inits(y, n, m, nclus, mu,grpids, nrandoms, (&maxval));
  }


  if(startmeth){
    if (reqsvd) {
      double **U, **V, *D, **xprime, *plrng, **theta, **ysamp;
      int j, d2, nm,ll;
      /*    FILE *fout;*/
      if(m < n) nm = m;
      else nm = n;

      MAKE_MATRIX(U, nm, m);
      MAKE_MATRIX(V, m, m);
      MAKE_VECTOR(D, m);
      MAKE_MATRIX(ysamp,nm,m);

 
      for(i = 0; i < nm; i++){
	for(j = 0; j < m; j++)	ysamp[i][j]=y[i][j];
      }

      d2 = svdd(ysamp, nm, m, D, U, V);
      
      FREE_MATRIX(ysamp);

      FREE_MATRIX(U);
      
      MAKE_MATRIX(xprime, n, (m - 1));

      for(i = 0; i < n; i++){
	for(j = 0; j < (m - 1); j++){
	  xprime[i][j]=0.;
	  for(ll = 0; ll < m; ll++) xprime[i][j] += y[i][ll] * V[ll][j];
	}
      }

      FREE_MATRIX(V);
      FREE_VECTOR(D);
      
      MAKE_MATRIX(theta, n, (m - 2));
      
      MAKE_VECTOR(plrng, 2);
      ind = get_polars_range_for_last_angl(n, m - 1, xprime, theta, plrng);
      FREE_MATRIX(xprime);
      
      if (ind) {
	for (i = 0; i < n; i++) {
	  if (theta[i][(m - 3)] >= plrng[1]) theta[i][(m - 3)] -= 2*PI;
	}
      }
      
      FREE_VECTOR(plrng);
      
      ind = starts_in_polar_domain(n, (m - 2), theta, nclus, ningrp, id);
      
      FREE_MATRIX(theta);

    }
    else {
      double *plrng, **theta;
      int ind;
      
      MAKE_MATRIX(theta, n, m - 1);
      
      MAKE_VECTOR(plrng, 2);
      ind = get_polars_range_for_last_angl(n, m, y, theta, plrng);
      
      if (ind) {
	for (i = 0; i < n; i++) {
	  if (theta[i][m - 2] >= plrng[1]) theta[i][m - 2] -= 2*PI;
	}
      }
      FREE_VECTOR(plrng);
      
      ind = starts_in_polar_domain(n, m - 1, theta, nclus, ningrp, id);
      
      FREE_MATRIX(theta);
      
    }
  }
  MAKE_MATRIX(ltsigma,nclus,m*(m+1)/2);

  temp = Inf;

  if(!ind){
    if(startmeth){
      i = initials(y, n, m, nclus, ningrp, mu, ltsigma, id);
      spherize(mu, nclus, m);
      temp =  obj_function(y, mu, id, n,m);
    }
    if(temp < maxval){
      maxval = temp;
      for(i=0;i<n;i++) grpids[i]=id[i];
    }
  }

  FREE_MATRIX(ltsigma);

  FREE_MATRIX(mu);
  FREE_VECTOR(id);

  return ind;
}

int kmeandirstarts(int n, int p, int nclus, double **x, double **Mu, 
		   int reqsvd, int startmeth, int nrandoms) 
{
  int i, ind = 0, *grpids, *ningrp;
  double **ltsigma;
  
  MAKE_VECTOR(ningrp, nclus);
  MAKE_VECTOR(grpids, n);
  MAKE_MATRIX(ltsigma,nclus,p*(p+1)/2);

  if(nclus > 1){
    ind = jumpstarters(n, p, x, nclus, ningrp, grpids, reqsvd, 
		       startmeth, nrandoms);
  }

  i = initials(x, n, p, nclus, ningrp, Mu, ltsigma, grpids);
  spherize(Mu, nclus, p);

  FREE_MATRIX(ltsigma);
  FREE_VECTOR(grpids);
  FREE_VECTOR(ningrp);
  return ind;
}

