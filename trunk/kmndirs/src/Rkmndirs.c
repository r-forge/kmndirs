#include <R.h>
#include <Rmath.h>

#include "array.h"

int run_kmndirs(double **X, int n, int p, double **mu, int K, int *id, 
		double **startingmu, int iter, double *normc, int *nc,
		int nrand);

void
R_kmndirs(double *x, int *nrx, int *ncx, int *k, int *nruns,
	  int *iterations, int *ids)
{
    double **data, *val;
    double **Mu, **startingmu, *normc;
    int i, j, *nc;

    /* Copy the R matrix data to C matrix. */
    data = (double **) R_alloc(*nrx, sizeof(double));
    for(i = 0; i < *nrx; i++) {
        data[i] = (double *) R_alloc(*ncx, sizeof(double));
        val = x + i;
        for(j = 0; j < *ncx; j++, val += *nrx)
            data[i][j] = *val;
    }

    MAKE_MATRIX(Mu, *k, *ncx);
    MAKE_MATRIX(startingmu, *k, *ncx);
    MAKE_VECTOR(normc, *k);
    MAKE_VECTOR(nc, *k);

    run_kmndirs(data, *nrx, *ncx, Mu, *k, ids,
		startingmu, *iterations, normc, nc, *nruns);

    FREE_MATRIX(Mu);
    FREE_MATRIX(startingmu);
    FREE_VECTOR(nc);
    FREE_VECTOR(normc);
}    
