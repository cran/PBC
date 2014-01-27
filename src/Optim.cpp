/*
 * Compute arg min of a function (C++) given
 *
 * Two methods:
 *
 * 1. Broyden–Fletcher–Goldfarb–Shanno algorithm (BFGS)
 * (based on 'bfgs' function of optim.c in R source code)
 *
 * 2. Limited-memory BFGS with bounds (L-BFGS-B)
 * (based on 'vmmin' function of optim.c in R source code)
 *
 */

#include "Optim.h"
#include <limits>

SEXP binMat;
double *ddata;
int nbObs;
int dimPara;
int dimData;
SEXP root;
SEXP f;
SEXP nIteration;
SEXP dxf;
SEXP dxdyf;
SEXP graf;
SEXP gradxf;
SEXP gradxdyf;
SEXP typef;
SEXP g;
SEXP dxg;
SEXP o;

void setParameters(SEXP m_binMat, SEXP m_x, SEXP m_root, SEXP m_f, SEXP m_nIteration, SEXP m_dxf, SEXP m_dxdyf,
                    SEXP m_graf, SEXP m_gradxf, SEXP m_gradxdyf, SEXP m_type, SEXP m_g, SEXP m_dxg, SEXP m_out, SEXP par){
  binMat = m_binMat;
  ddata = REAL(m_x);
  root = m_root;
  f = m_f;
  nIteration = m_nIteration;
  dxf = m_dxf;
  dxdyf = m_dxdyf;
  gradxf = m_gradxf;
  gradxdyf = m_gradxdyf;
  typef = m_type;
  g = m_g;
  dxg = m_dxg;
  o = m_out;
  dimPara = LENGTH(par);
  dimData = dimPara + 1;
  nbObs = (int) LENGTH(m_x) / dimData;
}

static double* fmin(int n, double *p)
{
  PBC pbc;
  NumericVector gradient(dimPara);
  double density = 0;
  NumericVector theta(n);
  for (int i = 0; i < n; i++)
    theta[i] = p[i];
  for (int i = 0; i < nbObs; i++){
	NumericVector tmpX(dimData);
	for (int j = 0; j < dimData; j++){
	  tmpX[j] = ddata[i * dimData + j];
    }
    pbc.parametersInit(binMat, theta, tmpX, root, f, nIteration, dxf,
                       dxdyf, graf, gradxf, gradxdyf, typef, g, dxg, o); // Parameters transmission
    NumericVector resultat;
    resultat = pbc.algo();
    for (int j = 0; j < dimPara; j++)
      gradient[j] += - resultat[j + 1] / resultat[0];
    density += - log(resultat[0]);
  }
  double* ret = new double[n + 1];
  ret[0] = density;
  for (int i = 1; i < n + 1; i ++) {
    ret[i] = gradient[i - 1];
  }
  return ret;
}

RcppExport SEXP lbfgsb(SEXP par, SEXP lower, SEXP upper,
                SEXP m_binMat, SEXP m_x, SEXP m_root, SEXP m_f, SEXP m_nIteration, SEXP m_dxf, SEXP m_dxdyf,
                SEXP m_graf, SEXP m_gradxf, SEXP m_gradxdyf, SEXP m_type, SEXP m_g, SEXP m_dxg, SEXP out){

	setParameters(m_binMat, m_x, m_root, m_f, m_nIteration, m_dxf,
                       m_dxdyf, m_graf, m_gradxf, m_gradxdyf, m_type, m_g, m_dxg, out, par);
	int npar = LENGTH(par); // Length of variable
	double* x = REAL(par); // Starting point
	double* l = REAL(lower); // Lower of variable
	double* u = REAL(upper); // Upper of variable
	double Fmin = 0.0; // Value minimum of function
	int lmm = 5;
	int trace = 0;
	int maxit = 100;
	double pgtol = 0.0;
	double factr = 1.0e+8;
	int nREPORT = 10;
	int* nbd = new int[npar];
	for (int i = 0; i < npar; i++)
	  nbd[i] = 2;
	int fail = 0;
	int fncount = 0;
	int grcount = 0;
	char msg[60];
    lbfgsb_modified(npar, lmm, x, l, u, nbd, &Fmin, fmin, &fail,
	    factr, pgtol, &fncount, &grcount,
	      maxit, msg, trace, nREPORT);
    NumericVector ret(npar);
    for (int i = 0; i < npar; i++)
      ret[i] = x[i];
	return ret;
}

RcppExport SEXP bfgs(SEXP par, SEXP m_binMat, SEXP m_x, SEXP m_root, SEXP m_f, SEXP m_nIteration, SEXP m_dxf, SEXP m_dxdyf,
                SEXP m_graf, SEXP m_gradxf, SEXP m_gradxdyf, SEXP m_type, SEXP m_g, SEXP m_dxg, SEXP out){

	setParameters(m_binMat, m_x, m_root, m_f, m_nIteration, m_dxf,
                       m_dxdyf, m_graf, m_gradxf, m_gradxdyf, m_type, m_g, m_dxg, out, par);
	int npar = LENGTH(par); // length of variable
	double* x = REAL(par); // starting point
	double Fmin = 0.0; // Value minimum of function
	int trace = 0;
	int maxit = 100;
	int fail = 0;
	int fncount = 0;
	int grcount = 0;
	int nREPORT = 10;
    int *mask;
    mask = (int *) R_alloc(npar, sizeof(int));
	for (int ii = 0; ii < npar; ii++) mask[ii] = 1;
	double abstol = numeric_limits<double>::max();
	double reltol = 1.490116e-08;
    vmmin_modified(npar, x, &Fmin, fmin, maxit, trace, mask,
      (-abstol), reltol, nREPORT, &fncount, &grcount, &fail);
    NumericVector ret(npar);
    for (int i = 0; i < npar; i++)
      ret[i] = x[i];
	return ret;
}
