#include <math.h>

void spherize(double **x, int nrow, int ncol)
{
  double xvar;
  int i, j;

  for (i=0; i<nrow; i++) {
    xvar = 0.0;
    for (j=0; j<ncol; j++) xvar += x[i][j] * x[i][j];
    if (xvar > 0) {
      for (j=0; j<ncol; j++) x[i][j] /= sqrt(xvar);
    }
  }

  return;
}
