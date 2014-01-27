/*
 *  C++ Header of derivatives  
 */

#ifndef DERIVATIVES_H
#define DERIVATIVES_H

#include <vector>
#include <iostream>

#include <math.h>
#include <Rcpp.h>

using namespace std;

SEXP compute(double x, double paraX, double y, double paraY, double theta, SEXP f, SEXP g);
double  Phi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type);
double  dxPhi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type);
double  dxdyPhi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type);
vector<double>  graPhi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type);
vector<double>  graDxPhi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type);
vector<double>  graDxDyPhi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type);

#endif // DERIAVATIVES_H
