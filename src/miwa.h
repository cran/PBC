/*
 * Include file for calculating orthant probabilities
 *
 * Stored in
 *   "orthant.h"
 *
 * Last modified:
 *   2003-05-06
 *
 * (c) 2001-2003 T. Miwa
 *
 */

/*
 * 
 * Modified for PBC package (R) (by Pham Van Trung)
 *
 */

#ifndef MIWA_H
#define MIWA_H

/* include R header files */
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
//#include <Rcpp.h>

#define MAXM    20    /* maximum dimension size  */
#define MAXGRD  4097  /* maximum number of grid points */

/* normal density function */
#define nrml_dn(X)  Rf_dnorm4(X, 0, 1, 0)

/* 1 - normal cumulative distribution function */
#define nrml_cd(X)  Rf_pnorm5(X, 0, 1, 1, 0)

/* Grid information
 * The values are calculated in gridcalc.c
 *   and passed to "orschm.c" throuth "orthant.c".
 */
struct GRID{
  int n;                /* even number of grid points */
  double z[MAXGRD];     /* grid points */
  double w[MAXGRD];     /* width, w[j]=z[j]-z[j-1] */
  double p[MAXGRD];     /* lower normal prob at z[j] */
  double d[MAXGRD];     /* normal density at z[j] */
  double w2[MAXGRD];    /* w squared */
  double w3[MAXGRD];    /* w cubed */
  double q[MAXGRD][4];  /* integral (x - z[j-1])^i * phi(x) */
};

void b_calc(int j, struct GRID *g, double *f, double *df,
                   double *b);
double dlt_f(int j, struct GRID *g,
                    double np, double nd, double dz,
                    double *b);
double orschm(int m, double *r, double *h, struct GRID *g);
double orthant(int m, double r[][MAXM][MAXM], double h[][MAXM],
               int *ncone, struct GRID *grid);
double nrml_lq(double p, double ueps, double peps, int *itr);
void gridcalc(struct GRID *g);
int checkall(int *vector,int length,int value);
double C_miwaTest(int steps, double* dcorr, double* dupper, double* dlower, int* infinvalue);
double calcul(double u1, double n1, double u2, double n2, double a);

#endif // MIWA_H
