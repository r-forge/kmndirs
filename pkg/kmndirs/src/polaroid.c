#include <stdio.h>
#include <stdlib.h>
#include "array.h"
#include <math.h>
#define PI 3.1415926535
#define NEG(x) ((x<0) ? 1 : 0)
#define SIGN(x) ((x) > 0 ? 1 : -NEG(x)) /* mimics the Splus sign function */


int polaroid(int p, double *x, double *R, double *theta)
/* This function calculates the polar coordinates of a p-variate vector x and 
   returns its norm, R and the (p-1)-dimensional angular coordinates theta. The
   function returns the index of the elements in x, after which all the 
   coordinates are zero. This means that the theta's from that index are 
   undefined.

   Written by xxxxxxxxxxxxx. */
{
	int i;
	double sum = 0.;
	
	for (i = 0; i < p; i++) sum += x[i]*x[i];
	
	(*R) = sqrt(sum);

	for (i = 0; i < (p - 1); i++) theta[i] = 0.;

	if (sum == 0) return 0;
	else {
		int k;
		double *y;

		MAKE_VECTOR(y, p - 1);    

		y[0] = sum;

		for (i = 1; i < (p - 1); i++) y[i] = y[i-1] - x[i-1]*x[i-1];

		for (k = 1; (k <= (p - 1)) && (y[k - 1] > 0); k++);

		for (i = 0; i < (k - 1); i++) {
			theta[i] =  acos( x[i] / sqrt(y[i]));
		}
		
		if (k == p) {
/*			if (SIGN(sin(theta[p - 2])) != SIGN(x[p - 1])) {*/
			if (NEG(x[p - 1])) {
/*				printf("\n\nhere we are %d %d %d\n\n", NEG(sin(theta[p - 2] / x[p - 1])), k, p);
				
				printf("x[%d] = %f, y[%d] = %f\n", p-1, x[p-1],  p-2, sqrt(y[p-2]));
				printf("\n cos = %e sin = %e x/r = %e\n",
				       cos(theta[p-2]), sin(theta[p - 2]),  x[p - 1] / sqrt(y[p - 2]));
*/
				theta[p - 2] = 2*PI - theta[p - 2];
/*				printf("\n cos = %e sin = %e x/r = %e\n", 
				       cos(theta[p-2]), sin(theta[p - 2]),  x[p - 1] / sqrt(y[p - 2]));
*/			}
		}

		FREE_VECTOR(y);

		return k;
	}
}

void euclid(int p, double *y, double R, double *theta)
/* This function reverses and calculates a p-variate x from its polar 
   coordinate representation. The input is R and the (p-1)-dimensional 
   angular coordinates theta. All angles  other than the last are assumed to 
   be between 0 and PI: the last (p-1)'th angle is between 0 and 2*PI.

   Written by xxxxxxxxxxxx */
{
	int i;
  
/* 
  for (i = 0; i < (p-1); i++) {
    y[i] = R * cos( theta[i] );
   printf("theta [%d] = %f y[%d] = %f\n", i, theta[i], i, y[i]);
    for (j = 0; j <= (i-1); j++) y[i] *= sin( theta[j] );
    printf("y[%d] = %f\n", i, y[i]);
  }
  y[p-1] = R;
  for (j = 0; j < (p-1); j++) y[p-1] *= sin( theta[j] ); */

  y[0] = R;

/*  for (i = 0; i < (p-1); i++) printf("%f ", theta[i]);*/

  for (i = 1; i < (p-1); i++) y[i] = y[i - 1] * sin( theta[ i - 1 ] );
  y[p - 1] = y[ p - 2 ] * sin( theta[ p - 2 ]);
  for (i = 0; i < (p-1); i++) y[i] *= cos( theta[ i ]);

  return;
}
