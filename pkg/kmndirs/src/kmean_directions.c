#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define BIG 1e+40
#include "array.h"
/* Code based on Applied Statistics algorithms (C) Royal Statistical Society
   1979. Adapted for C by xxxxxxxxxxxxxxxxx */

  
/* identical to kmns.c except for indx which is a number here */


double unitnorm_sphere(double *m,int lengthm);

void dir_optra(double **a,int m,int n,double **c,int k,int *ic1,int *ic2,
	       int *nc, double *an1,double *an2,int *ncp,double *d,int *itran,
	       int *live, int *indx,double *normc);

void dir_qtran(double **a,int m,int n,double **c,int k,int *ic1,int *ic2,
	       int *nc,double *an1,double *an2,int *ncp,double *d,int *itran,
	       int *indx,double *normc);

void kmean_dirs(double **a, int m, int n, double **c, int k, int *ic1,int *nc,
		int iter,double *wss,int *ifault,double *normc)
{
  int i,j,l,ij,il,indx,*ic2,*ncp,*itran,*live;
  double temp,da,db,dc,dt[2],*an1,*an2,*d;
  
  /* ALGORITHM AS 136  APPL. STATIST. (1979) VOL.28, NO.1 
     Divide M points in N-dimensional space into K clusters so that the within 
     cluster sum of squares is minimized. */ 
  
  MAKE_VECTOR(ic2,m);

  *ifault = 3;
  if ((k <=1) || (k >= m)) {  /* check for number of clusters */
    return;
  }
  *ifault = 0;

  /* For each point i, find its two closest centres, ic1(i) and ic2(i).
     Assign it to ic1(i). */

  for (i=0;i<m;i++) {
    ic1[i] = 0;
    ic2[i] = 1;
    for (il=0;il<2;++il) {
      dt[il]=0.;
      for (j=0;j<n;++j) {
	da=a[i][j]-c[il][j];
	dt[il]+=da*da;
      }
    }
    if (dt[0] >= dt[1]) {
      ic1[i] = 1;
      ic2[i] = 0;
      temp=dt[0];
      dt[0] = dt[1];
      dt[1] = temp;
    }
    for (l=2;l<k;++l) {
      db=0.;
      for (j=0;j<n;++j) {
	dc=a[i][j]-c[l][j];
	db+=dc*dc;
	if (db>dt[1]) {
	  j=n; /*go to the end of the loop -- end of story*/
	}
      }
      if (db<dt[1]) {
	if (db>=dt[0]) {
	  dt[1] = db;
	  ic2[i] = l;
	}
	else {
	  dt[1] = dt[0];
	  ic2[i] = ic1[i];
	  dt[0] = db;
	  ic1[i] = l;
	}
      }
    }
  }

/* Update cluster centers to be the average of points contained within them. */

  MAKE_VECTOR(an1,k);
  MAKE_VECTOR(an2,k);
  MAKE_VECTOR(d,m);
  MAKE_VECTOR(itran,k);
  MAKE_VECTOR(live,k);
  MAKE_VECTOR(ncp,k);

  for (l=0;l<k;++l) {
    nc[l] = 0;
    for (j=0;j<n;++j) {
      c[l][j]=0.;
    }
  }
  for (i=0;i<m;++i) {
    nc[ic1[i]]++;
    for (j=0;j<n;++j) {
      c[ic1[i]][j]+=a[i][j];
    }
  }

  /*  Check to see if there is any empty cluster at this stage */

  for (l=0;l<k;++l) {
    if (nc[l] == 0) {
      *ifault=1;
      return;
    }
    for (j=0;j<n;++j) {
      c[l][j]/=nc[l];
    }

    normc[l]=unitnorm_sphere(c[l],n);

/*    printf("l=%d  %lf\n",l,normc[l]);*/

    /* Initialize an1, an2, itran & ncp: 
       an1[l]=nc[l]/(nc[l]-1); an2[l]=nc[L]/(nc[l]+1)
       itran[l] = 1 if cluster l is updated in the quick-transfer stage, 0 ow.
       In the optimal-transfer stage, ncp(l) stores the step at which cluster
       l is last updated. In the quick-transfer stage, ncp(l) stores the step 
       at which cluster l is last updated plus m. */
    
    /*    an2[l] =nc[l]/(nc[l]+1.0);*/
    an2[l] = 1.*nc[l]+1. - 1.*nc[l]*nc[l]*normc[l]*normc[l]/(nc[l]+1.) 
      - 1./(nc[l]+1.);  
    an1[l] = BIG;


    if (nc[l]>1) {
      /*      an1[l] =nc[l]/(nc[l]-1.0);*/
      an1[l] = 1.*nc[l]-1. - 1.*nc[l]*nc[l]*normc[l]*normc[l]/(nc[l]-1.) 
	- 1./(nc[l]-1.);
    }
    itran[l] = 1;
    ncp[l] = -1;
  }

  indx=0;
  for (ij=0;ij<iter;++ij) {
    /* In this stage, there is only one pass through the data. Each point is 
       re-allocated, if necessary, to the cluster that will induce the maximum
       reduction in within-cluster sum of squares. */

    dir_optra(a,m,n,c,k,ic1,ic2,nc,an1,an2,ncp,d,itran,live,&indx,normc);
  
    /* Stop if no transfer took place in the last m optimal transfer steps. */
    if (indx == m) {
      ij=iter;
    }
    else { /* Each point is tested in turn to see if it should be re-allocated 
	      to the cluster to which it is most likely to be transferred, 
	      ic2[i], from its present cluster, ic1[i]. Loop through the data 
	      until no further change is to take place. */

      dir_qtran(a,m,n,c,k,ic1,ic2,nc,an1,an2,ncp,d,itran,&indx,normc);

      if (k==2) { /* k=2: no need to re-enter the optimal transfer stage*/
	ij=iter;
      }
      else {	/* ncp has to be set to 0 before entering optra. */
	for (l=0;l<k; ++l) {
	  ncp[l] = 0;
	}
      }
    }
  } 
  if ((indx!=m) && (k!=2)) { 
    *ifault = 2; /* iterations exceeded: may indicate unforeseen looping */
  }

  /* Compute within-cluster SS for each cluster. Not needed?*/
  /*
  for (l=0;l<k;++l) {
    wss[l]=0.;
    for (j=0;j<n;++j) {
      c[l][j]=0.;
    }
  }
  for (i=0;i<m;++i) {
    for (j=0;j<n;++j) {
      c[ic1[i]][j]+=a[i][j];
    }
  }
  for (j=0;j<n;++j) {
    for (l=0;l<k;++l) {
      c[l][j] /= (double)nc[l];
    }
  }

  for(i=0;i<k;i++){
    unitize(c[i],n);
  }
  
  for(j=0;j<n;++j){
    for (i=0;i<m;++i) {
      da = a[i][j]-c[ic1[i]][j];
      wss[ic1[i]]+=da*da;
    }
  }
*/
  FREE_VECTOR(ic2);
  FREE_VECTOR(an1);
  FREE_VECTOR(d);
  FREE_VECTOR(itran);
  FREE_VECTOR(live);
  FREE_VECTOR(ncp);
  FREE_VECTOR(an2);
  return;
} 

void dir_optra(double **a,int m,int n,double **c,int k,int *ic1,int *ic2,int *nc, 
	   double *an1,double *an2,int *ncp,double *d,int *itran,int *live,
	       int *indx, double *normc)
{
  int i,j,l,l1,l2,ll,flag;
  /*  double r2,da,db,dc,dd,de,df,rr,al1,al2,alt,alw;*/
  double r2,da,dc,de;
 /* ,rr,*oldl1,*oldl2,tmp1,tmp2*/
  /* ALGORITHM AS 136.1  APPL. STATIST. (1979) VOL.28, NO.1 
     This is the optimal transfer stage. Each point is re-allocated, if 
     necessary, to the cluster that will induce a maximum reduction in the 
     within-cluster sum of squares. 
     If cluster L is updated in the last quick-transfer stage, it belongs to 
     the live set throughout this stage.   Otherwise, at each step, it is not 
     in the live set if it has not been updated in the last M optimal transfer
     steps. */

  for (l=0;l<k;++l) {
    if (itran[l] == 1) {
      live[l]=m+1;
    }
  }
  for (i=0;i<m;++i) {
    (*indx)++;
    l1 = ic1[i];
    l2 = ic2[i];
    ll = l2;
    /* If point I is the only member of cluster L1, no transfer. */
    if (nc[l1] != 1) {
      /* If L1 has not yet been updated, no need to re-compute D(I).*/

      if (ncp[l1] != 0) {
	de=0.;
	/*	for (j=0;j<n;++j) {
	  df=a[i][j]-c[l1][j];
	  de+=df*df;
	}
	d[i]=de*an1[l1];*/
	for(j=0;j<n;++j){
	  de+= a[i][j]*c[l1][j];
	}
	d[i] = an1[l1]+2.*nc[l1]*normc[l1]*de/(nc[l1]-1.);
/*printf("%lf\n",normc[l1]);*/
      }
       /* Find the cluster with minimum R2. */
      da=0.;
      /*
	for (j=0;j<n;++j) {
	db=a[i][j]-c[l2][j];
	da+=db*db;
	}
	r2=da*an2[l2];
      */
      for(j=0;j<n;++j){
	da += a[i][j]*c[l2][j];
      }
      r2 = an2[l2]-2.*nc[l2]*normc[l2]*da/(nc[l2]+1.);
/*printf("%lf\n",r2);*/
      for (l=0;l<k;++l) {
	/* If I >= LIVE(L1), then L1 is not in the live set.   If this is 
	   true, we only need to consider clusters that are in the live set 
	   for possible transfer of point I.   Otherwise, we need to consider
	   all possible clusters. */
	if (((i>=(live[l1]-1)) && (i>=(live[l]-1))) || (l == l1) || (l == ll)) {
	}
	else {
	  /*	  rr=r2/an2[l];*/
/*	  rr = r2-an2[l];*/
	  dc=0.;
	  flag=0;
	  /*	  j=0;*/

	  for(j=0;j<n;j++){
	    dc+=a[i][j]-c[l][j];
	  }
	  dc*= -2.*nc[l]*normc[l]/(nc[l]+1.);
	  dc += an2[l];
	  if(dc>= r2) {
	    flag=1;
	  }

	  /*	  while ((flag==0) && (j<n)) {
		  dd=a[i][j]-c[l][j];
		  dc+=dd*dd;
		  if (dc>=rr) {
		  flag=1;
		  }
		  j++;
		  }*/

	  if (flag==0) {
	    /*	     r2=dc*an2[l];*/
/*	    r2 = dc + an2[l];*/
            r2=dc;
	    l2=l;
	  }
	}
      }
      if (r2>=d[i]) {/* If no transfer is necessary, L2 is the new IC2(I). */
	ic2[i] = l2;
      }
      else {/* Update cluster centres, LIVE, NCP, AN1 & AN2 for clusters L1 and
	       L2, and update IC1(I) & IC2(I). */
	(*indx)=0;
	live[l1]=m+i+1;
	live[l2]=m+i+1;
	ncp[l1]=i+1;
	ncp[l2]=i+1;

/*	MAKE_VECTOR(oldl1,n);
	MAKE_VECTOR(oldl2,n);*/

	/* Update means  */
/*	for(j=0;j<n;j++){
	  oldl1[j]=c[l1][j];
	  oldl2[j]=c[l2][j];
	}*/



	for(j=0;j<n;j++){
	  c[l1][j]= (double)nc[l1]*normc[l1]*c[l1][j]-a[i][j];
	  c[l2][j]=(double)nc[l2]*normc[l2]*c[l2][j]+a[i][j];
	  c[l1][j]/=(1.*nc[l1]-1.);
	  c[l2][j]/=(1.*nc[l2]+1.);
	}


	/*  Update norms */
/*	tmp1=0.;
	tmp2=0.;
	for(j=0;j<n;j++){
	  tmp1+=a[i][j]*oldl1[j];
	  tmp2+=a[i][j]*oldl2[j];
	}


	normc[l1] = 1.*nc[l1]*nc[l1]*normc[l1]*normc[l1] - nc[l1]*normc[l1]*tmp1 + 1.;
	normc[l1]=sqrt(normc[l1])/(1.*nc[l1]-1.);
	normc[l2] = 1.*nc[l2]*nc[l2]*normc[l2]*normc[l2] + nc[l2]*normc[l2]*tmp2 + 1.;
	normc[l2]=sqrt(normc[l2])/(1.*nc[l1]+1.);

	for(j=0;j<n;j++){
	  c[l1][j]/=normc[l1];
	  c[l2][j]/=normc[l2];
	}
*/



/*	FREE_VECTOR(oldl1);
	FREE_VECTOR(oldl2);*/


	nc[l1]--;
	nc[l2]++;


/*	an2[l2] = 1.*nc[l2]+1. - 1.*nc[l2]*nc[l2]*normc[l2]*normc[l2]/(nc[l2]+1.) 
	  - 1./(nc[l2]+1.);  
	an2[l1] = 1.*nc[l1]+1. - 1.*nc[l1]*nc[l1]*normc[l1]*normc[l1]/(nc[l1]+1.) 
	  - 1./(nc[l1]+1.);
	an1[l2] = 1.*nc[l2]-1. - 1.*nc[l2]*nc[l2]*normc[l2]*normc[l2]/(nc[l2]-1.) 
	  - 1./(nc[l2]-1.);
	an1[l1] = 1.*nc[l1]-1. - 1.*nc[l1]*nc[l1]*normc[l1]*normc[l1]/(nc[l1]-1.) 
	  - 1./(nc[l1]-1.);
*/
	ic1[i]=l2;
	ic2[i]=l1;


normc[l1]=unitnorm_sphere(c[l1],n);
normc[l2]=unitnorm_sphere(c[l2],n);

	an2[l2] = 1.*nc[l2]+1. - 1.*nc[l2]*nc[l2]*normc[l2]*normc[l2]/(nc[l2]+1.) 
	  - 1./(nc[l2]+1.);  
	an2[l1] = 1.*nc[l1]+1. - 1.*nc[l1]*nc[l1]*normc[l1]*normc[l1]/(nc[l1]+1.) 
	  - 1./(nc[l1]+1.);
	an1[l2] = 1.*nc[l2]-1. - 1.*nc[l2]*nc[l2]*normc[l2]*normc[l2]/(nc[l2]-1.) 
	  - 1./(nc[l2]-1.);
	an1[l1] = 1.*nc[l1]-1. - 1.*nc[l1]*nc[l1]*normc[l1]*normc[l1]/(nc[l1]-1.) 
	  - 1./(nc[l1]-1.);

/*printf("%lf\n",an1[l2]);*/
      }
    }
    if ((*indx)==m) {
      return;
    }
  }

  for (l=0;l<k;++l) { /* ITRAN(L) = 0 before entering QTRAN. Also, LIVE(L) has
			 to be decreased by M before re-entering OPTRA. */ 
    itran[l]=0;
    live[l]-=m;
  }
  
  
  return;
}

void dir_qtran(double **a,int m,int n,double **c,int k,int *ic1,int *ic2,
	       int *nc,double *an1,double *an2,int *ncp,double *d,int *itran,
	       int *indx,double *normc)
{
  int i,j,l1,l2,icoun,istep,flag,iflag;
  double r2,da,dd;
/*,*oldl1,*oldl2,tmp1,tmp2;*/
  
  /* ALGORITHM AS 136.2  APPL. STATIST. (1979) VOL.28, NO.1 
   This is the quick transfer stage. IC1(I) is the cluster which point I 
   belongs to. IC2(I) is the cluster which point I is most likely to be 
   transferred to. For each point I, IC1(I) & IC2(I) are switched, if 
   necessary, to reduce within-cluster sum of squares.  The cluster centres 
   are updated after each step.
   In the optimal transfer stage, NCP(L) indicates the step at which 
   cluster L is last updated.   In the quick transfer stage, NCP(L) is 
   equal to the step at which cluster L is last updated plus M. */
  icoun=0;
  istep=0;
  iflag=0;
/*  while (iflag==0) {*/
  for (i=0;i<m;i++) {
    ++icoun;
    ++istep;
    l1 = ic1[i];
    l2 = ic2[i];
    
    /*  If point I is the only member of cluster L1, no transfer. */
    if (nc[l1] != 1) {/* If ISTEP > NCP(L1), no need to re-compute distance 
			 from point I to cluster L1.   Note that if cluster
			 L1 is last updated exactly M steps ago, we still
			 need to compute the distance from point I to 
			 cluster L1. */
      if (istep<=ncp[l1]) {
	da=0.;
	/*	for (j=0;j<n;++j) {
	  db = a[i][j]-c[l1][j];
	  da += db * db;
	}
	d[i]=da*an1[l1];*/

	for(j=0;j<n;++j){
	  da+= a[i][j]*c[l1][j];
	}
	d[i] = an1[l1]+2.*nc[l1]*normc[l1]*da/(nc[l1]-1.);

	/* If ISTEP >= both NCP(L1) & NCP(L2) there will be no transfer 
	   of point I at this step. */
      }
      if ((istep<ncp[l1]) || (istep < ncp[l2])) {
	r2 = d[i] - an2[l2];
	/*	r2=d[i]/an2[l2];*/
	dd=0.;
	flag=0;
	/*	j=0;
		while ((j<n) && (flag==0)) {
		de=a[i][j]-c[l2][j];
		dd+=de*de;
		if(dd>= r2) {
		flag=1;
		}
		j++;
		}*/

	for(j=0;j<n;j++){
	  dd+=a[i][j]-c[l2][j];
	}
	dd*= -2.*nc[l2]*normc[l2]/(nc[l2]+1.);

	if(dd>= r2) {
	  flag=1;
	}

	if (flag==0) {
	  /* Update cluster centres, NCP, NC, ITRAN, AN1 & AN2 for clusters L1
	     & L2. Also update IC1(I) & IC2(I). Note that if any updating 
	     occurs in this stage, INDX is set back to 0. */ 
	  icoun = 0;
	  (*indx) = 0;
	  itran[l1]=1;
	  itran[l2]=1;
	  ncp[l1]=istep+m;
	  ncp[l2]=istep+m;

/*	  MAKE_VECTOR(oldl1,n);
	  MAKE_VECTOR(oldl2,n);
*/	  
	  /* Update means  */
/*	  for(j=0;j<n;++j){
	    oldl1[j]=c[l1][j];
	    oldl2[j]=c[l2][j];
	  }*/

	  for(j=0;j<n;++j){
	    c[l1][j]=nc[l1]*normc[l1]*c[l1][j]-a[i][j];
	    c[l2][j]=nc[l2]*normc[l2]*c[l2][j]+a[i][j];
  	    c[l1][j]/=(1.*nc[l1]-1.);
            c[l2][j]/=(1.*nc[l2]+1.);
	  }
	  /*  Update norms */
/*	  tmp1=0.;
	  tmp2=0.;
	  for(j=0;j<n;j++){
	    tmp1+=a[i][j]*oldl1[j];
	    tmp2+=a[i][j]*oldl2[j];
	  }
	  normc[l1] = 1.*nc[l1]*nc[l1]*normc[l1]*normc[l1] - nc[l1]*normc[l1]*tmp1 + 1.;
	  normc[l1]=sqrt(normc[l1]);
	  normc[l2] = 1.*nc[l2]*nc[l2]*normc[l2]*normc[l2] + nc[l2]*normc[l2]*tmp2 + 1.;
	  normc[l2]=sqrt(normc[l2]);
	  
	  for(j=0;j<n;j++){
	    c[l1][j]/=normc[l1];
	    c[l2][j]/=normc[l2];
	  }
	  
	  FREE_VECTOR(oldl1);
	  FREE_VECTOR(oldl2);*/
	  
  
	  nc[l1]--;
	  nc[l2]++;

normc[l1]=unitnorm_sphere(c[l1],n);
normc[l2]=unitnorm_sphere(c[l2],n);

	  an2[l2] = 1.*nc[l2]+1. - 1.*nc[l2]*nc[l2]*normc[l2]*normc[l2]/(nc[l2]+1.) 
	    - 1./(nc[l2]+1.);  
	  an2[l1] = 1.*nc[l1]+1. - 1.*nc[l1]*nc[l1]*normc[l1]*normc[l1]/(nc[l1]+1.) 
	    - 1./(nc[l1]+1.);
	  an1[l2] = 1.*nc[l2]-1. - 1.*nc[l2]*nc[l2]*normc[l2]*normc[l2]/(nc[l2]-1.) 
	    - 1./(nc[l2]-1.);
	  an1[l1] = 1.*nc[l1]-1. - 1.*nc[l1]*nc[l1]*normc[l1]*normc[l1]/(nc[l1]-1.) 
	    - 1./(nc[l1]-1.);
	  ic1[i] = l2;
	  ic2[i] = l1;
	}
      }
    }
    if (icoun==m) return;
  }  
/*  }*/
}

/* This function is used to force a vector to be a unit vector with mean zero.
i.e sum(m)=0 and sum(m^2)=1 
  IMPORTANT! THIS CHANGES THE VALUES IN THE INPUT VECTOR. 
*/

double unitnorm_sphere(double *m,int lengthm)
{

  int i;
  double sum;
  sum=0.0;

  for(i=0;i<lengthm;i++){
    sum+=m[i]*m[i];
  }

  for(i=0;i<lengthm;i++){
    m[i]/=sqrt(sum);
  }
  
  return sqrt(sum);
}


/* This function standardizes the rows of x so that they have zero mean 
   and lie on the unit ncol-dimensional sphere */

/*
void standardize(double **x, int nrow, int ncol)
{
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
*/


