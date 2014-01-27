#ifndef OPTIM_H
#define OPTIM_H

#include <Rcpp.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "PBC.h"
#include "Derivatives.h"
#include "lbfgsb.h"

typedef double optimfn(int, double *);
typedef double* optimgr(int, double *);
typedef double* optim(int, double *);

#endif //OPTIM_H
