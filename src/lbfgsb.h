#ifndef LBFGSB_H
#define LBFGSB_H

#include <R.h>
#include <math.h>
#include <float.h> /* for DBL_EPSILON */
#include <string.h>
#include <R_ext/RS.h> /* for F77_CALL */
#include <R_ext/Linpack.h>
#include <R_ext/Print.h> /* Rprintf */

#include "Optim.h"

typedef double optimfn(int, double *);
typedef double* optimgr(int, double *);
typedef double* optim(int, double *);

void lbfgsb_modified(int n, int m, double *x, double *l, double *u, int *nbd,
	    double *Fmin, optim fmin, int *fail,
	    double factr, double pgtol,
	    int *fncount, int *grcount, int maxit, char *msg,
	    int trace, int nREPORT);
void vmmin_modified(int n0, double *b, double *Fmin, optim fmin,
      int maxit, int trace, int *mask,
      double abstol, double reltol, int nREPORT,
      int *fncount, int *grcount, int *fail);

#endif // LBFGSB_H
