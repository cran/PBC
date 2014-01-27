/*
 *  Derivatives.cpp : Get derivatives from R for computing PBC
 */

#include <iostream>

#include <math.h>
#include <Rmath.h>

#include "Derivatives.h"
#include "miwa.h"
#include "Gradient.h"

using namespace std;
using namespace Rcpp;

const int n = 500;

// Set parameters for derivative functions in R 
SEXP compute(double x, double paraX, double y, double paraY, double a, SEXP f, SEXP g){
  Function gFunc("margin");
  NumericVector PARAX = wrap(paraX);
  NumericVector PARAY = wrap(paraY);
  NumericVector X = gFunc(PARAX, wrap(x), g); // Convert from double in C++ to SEXP in Rcpp
  NumericVector Y = gFunc(PARAY, wrap(y), g);
  NumericVector A = wrap(a);
  Function func("compute"); // Call function "compute" from R
  NumericVector ret = func(X, Y, A, f); // Return value of function 
  return ret;
}

// Value of function
double Phi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type){
  double result = 0;
  double a = theta[0];
  if (type == 0){ // Distribution given
    NumericVector ret = compute(x, paraX, y, paraY, a, f, g);
    result = ret[0];
  }
  if (type == 1){// Gumbel distribution
    result = exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a));
  }
  if (type == 2){ // Gumbel copula
    a = (double) 1/a;
    result = exp(-pow(pow(-log(pow(x, 1/paraX)), 1/a) + pow(-log(pow(y, 1/paraY)), 1/a), a));
  }
  if (type == 3){ // Farlie-Gumbel-Morgenstern copula
    result = pow(x, 1/paraX) * pow(y, 1/paraY) * (1 + a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY)));
  }
  if (type == 4){ // Frank copula
    result = - log(1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1) / (exp(-a) - 1)) / a;
  }
  if (type == 5) { // Gaussian copula
    /*Function func("phi.norm");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A = wrap(a);
    NumericVector ret = func(X, PARAX, Y, PARAY, A); // Call function "phi.norm" from R
    result = ret[0];
    * */    
    result = calcul(x, paraX, y, paraY, a); // Call function calcul from miwa.cpp
  }
  if (type == 6) { // Ali-Mikhail-Haq copula
    result = pow(x, 1/paraX) * pow(y, 1/paraY) / (1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY)));
  }
  if (type == 7) { // Joe copula
    result = 1 - pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a)
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a);
  }
  if (type == 8){ // Student copula
    Function func("phi.student");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A(2);
    A[0] = theta[0];
    A[1] = theta[1];
    NumericVector ret = func(X, PARAX, Y, PARAY, A); // Call function "phi.student" from R
    result = ret[0];    
  }
  return result;
}

// Partial derivative of one variable
double  dxPhi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type){
  double result = 0;
  double a = theta[0];
  if (type == 0){ // Distribution given
    NumericVector ret = compute(x, paraX, y, paraY, a, f, g);
    Function dxgFunc("margin");
    NumericVector PARAX = wrap(paraX);
    NumericVector X = wrap(x);
    NumericVector deri = dxgFunc(PARAX, X, dxg);
    result = deri[0] * ret[0];
  }
  if (type == 1){ // Gumbel distribution
    result = exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) 
         * pow(x, -1/a - 1) * pow(paraX, -1/a) 
         * pow(pow(paraY * y, -1/a) + pow(paraX * x, -1/a), a - 1); 
  }
  if (type == 2){ // Gumbel copula
    a = (double) 1/a;
    result = - pow(-log(x) / paraX, 1/a) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a))
             * pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a - 1) / (x * log(x));
  }
  if (type == 3){ // Farlie-Gumbel-Morgenstern copula
    result = pow(x, 1/paraX - 1) * pow(y, 1/paraY) * (1 + a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))) / paraX
           - a * pow(x, 2/paraX - 1) * pow(y, 1/paraY) * (1 - pow(y, 1/paraY)) / paraX;
  }
  if (type == 4){ // Frank copula
    result = pow(x, 1/paraX - 1) * exp(-a * pow(x, 1/paraX)) * (exp(-a * pow(y, 1/paraY)) - 1) 
           / ((exp(-a) - 1) * paraX * (1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1)));
  }
 if (type == 5) { // Gaussian copula
    /*Function func("dxphi.norm");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A = wrap(a);
    NumericVector ret = func(X, PARAX, Y, PARAY, A); // Call function "dxphi.norm" from R
    result = ret[0];
    * */
    double xx[n+1];
    double eps;
    eps = (double) pow(y,(double)1/paraY)/n;
    xx[0] = 0;
    for(int i = 1; i < n+1; i++) {
      xx[i] = xx[i-1] + eps;
    }
    // Compute the integral
    double summand = 0;
    for(int i=0; i<n; i++){
      double mean = (xx[i+1] + xx[i]) / 2;
      double val1 = Rf_qnorm5(mean, 0.0, 1.0, 1, 0);
      double val2 = Rf_qnorm5(pow(x, (double) 1/paraX), 0.0, 1.0, 1, 0);
      double resultat = exp((pow(a * val1, 2) + pow(a * val2, 2) - 2 * a * val1 * val2)
                     / (-2 + 2 * a * a)) / pow((1 - a * a),0.5);
      summand += resultat * eps;
    }
    summand *= pow(x, (double)1/paraX - 1) / paraX;
    result = summand;    
  }
  if (type == 6) { // Ali-Mikhail-Haq copula
    result = pow(x, 1/paraX - 1) * pow(y, 1/paraY) / (paraX * (1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))))
           - a * pow(x, 2/paraX - 1) * pow(y, 1/paraY) * (1 - pow(y, 1/paraY)) 
           / (paraX * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 2));
  }
  if (type == 7) { // Joe copula
    result = - pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * (pow(1 - pow(y, 1/paraY), a) - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 1) / paraX;
  }
  if (type == 8){ // Student copula
    Function func("dxphi.student");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A(2);
    A[0] = theta[0];
    A[1] = theta[1];
    NumericVector ret = func(X, PARAX, Y, PARAY, A); // Call function "dxphi.student" from R
    result = ret[0];
  }
  return result;
}

// Partial derivative of two variables
double  dxdyPhi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type){
  double result = 0;
  double a = theta[0];
  if (type == 0){ // Distribution given
    NumericVector ret = compute(x, paraX, y, paraY, a, f, g);
    Function dxgFunc("margin");
    NumericVector PARAX = wrap(paraX);
    NumericVector X = wrap(x);
    NumericVector deri = dxgFunc(PARAX, X, dxg);
    NumericVector PARAY = wrap(paraY);
    NumericVector Y = wrap(y);
    NumericVector deri1 = dxgFunc(PARAY, Y, dxg);
    result = deri[0] * deri1[0] * ret[0];
  }
  if (type == 1){ // Gumbel distribution
    result = exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a))  
         * pow(pow(paraY * y, -1/a) + pow(paraX * x, -1/a), 2 * a - 2)
         * pow(x, -1/a - 1) * pow(paraX, -1/a)
         * pow(y, -1/a - 1) * pow(paraY, -1/a) 
         - (1 - 1/a) * pow(pow(paraY * y, -1/a) + pow(paraX * x, -1/a), a - 2)
         * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) 
         * pow(x, -1/a - 1) * pow(paraX, -1/a) 
         * pow(y, -1/a - 1) * pow(paraY, -1/a);
  }
  if (type == 2){ // Gumbel copula
    a = (double) 1/a;
    result =  pow(-log(x) / paraX, 1/a) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a))
           *  pow(-log(y) / paraY, 1/a) * pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), 2 * a - 2)
           / (x * log(x) * y * log(y))
           - (a - 1) * pow(-log(x) / paraX, 1/a) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a))
           *  pow(-log(y) / paraY, 1/a) * pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a - 2)
           / (a * x * log(x) * y * log(y));
  }
  if (type == 3){ // Farlie-Gumbel-Morgenstern copula
    result = a * pow(x, 2/paraX - 1) * pow(y, 2/paraY - 1) / (paraX * paraY)
           - a * pow(x, 1/paraX - 1) * (1 - pow(x, 1/paraX)) * pow(y, 2/paraY - 1) / (paraX * paraY)
           - a * pow(y, 1/paraY - 1) * (1 - pow(y, 1/paraY)) * pow(x, 2/paraX - 1) / (paraX * paraY)
           + pow(x, 1/paraX - 1) * pow(y, 1/paraY -  1) * (1 + a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))) / (paraX * paraY);
  }
  if (type == 4){ // Frank copula
    result = a * pow(x, 1/paraX - 1) * pow(y, 1/paraY - 1) * exp(-a * pow(x, 1/paraX) - a * pow(y, 1/paraY)) * (exp(-a * pow(y, 1/paraY)) - 1) * (exp(-a * pow(x, 1/paraX)) - 1) 
           / (pow((exp(-a) - 1), 2) * paraX * paraY * pow((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1)), 2))
           - a * pow(x, 1/paraX - 1) * pow(y, 1/paraY - 1) * exp(-a * pow(x, 1/paraX) - a * pow(y, 1/paraY)) 
           / ((exp(-a) - 1) * paraX * paraY * (1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1)));
  }
  if (type == 5) { // Gaussian copula
    /*Function func("dxdyphi.norm");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A = wrap(a);
    NumericVector ret = func(X, PARAX, Y, PARAY, A); // Call function "dxdyphi.norm" from R
    result = ret[0];
    * */
    double val1 = Rf_qnorm5(pow(x, 1/paraX), 0.0, 1.0, 1, 0);
    double val2 = Rf_qnorm5(pow(y, 1/paraY), 0.0, 1.0, 1, 0);
    result = pow(x, 1/paraX - 1) * pow (y, 1/paraY - 1) / (paraX * paraY);
    result *= exp((pow(a * val1, 2) + pow(a * val2, 2) - 2 * a * val1 * val2)
           / (-2 + 2 * a * a)) / pow((1 - a * a), 0.5);          
  }
  if (type == 6) { // Ali-Mikhail-Haq copula
    result = - a * pow(x, 2/paraX - 1) * pow(y, 1/paraY - 1) * (1 - pow(y, 1/paraY)) 
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 2))
           + 2 * a * a * pow(x, 2/paraX - 1) * pow(y, 2/paraY - 1) * (1 - pow(y, 1/paraY)) * (1 - pow(x, 1/paraX)) 
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 3))
           + pow(x, 1/paraX - 1) * pow(y, 1/paraY - 1)
           / (paraX * paraY * (1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))))
           + a * pow(x, 2/paraX - 1) * pow(y, 2/paraY - 1)
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 2))
           - a * pow(x, 1/paraX - 1) * pow(y, 2/paraY - 1) * (1 - pow(x, 1/paraX)) 
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 2));
  }
  if (type == 7) { // Joe copula
    result = a * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 1) / (paraX * paraY)
           - (1 - a) * pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1) * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1)
           * (pow(1 - pow(x, 1/paraX), a) - 1) * (pow(1 - pow(y, 1/paraY), a) - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 2) / (paraX * paraY);
  }
  if (type == 8){ // Student copula
    Function func("dxdyphi.student");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A(2);
    A[0] = theta[0];
    A[1] = theta[1];
    NumericVector ret = func(X, PARAX, Y, PARAY, A); // Call function "dxdyphi.student" from R
    result = ret[0];
  }
  return result;
}

// Gradient of a
vector<double>  graPhi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type){
  double result;
  vector<double> res;
  double a = theta[0];
  if (type == 0){ // Distribution given
    NumericVector ret = compute(x, paraX, y, paraY, a, f, g);
    result = ret[0];
    res.push_back(result);
  }
  if (type == 1){ // Gumbel distribution
    result =  - pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)
         * (log(pow(paraX * x, -1/a) + pow(paraY * y, -1/a)) 
         + a * (log(paraY * y) * pow(paraY * y, -1/a) * pow(a,-2) + log(paraX * x) * pow(paraX * x, -1/a) * pow(a,-2)) 
         / (pow(paraX * x, -1/a) + pow(paraY * y, -1/a)))
         * exp(-(pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)));
    res.push_back(result);  
  }
  if (type == 2){ // Gumbel copula
    double a_tmp = a;
    a = (double) 1/a;
    result = - exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a)) * pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a)
           *   (log(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)) - ((pow(-log(y)/paraY, 1/a)*log(-log(y)/paraY)) + (pow(-log(x)/paraX, 1/a)*log(-log(x)/paraX)))
           /   (a * (pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a))));
    result = result * (-1/(a_tmp * a_tmp));
    res.push_back(result);
  }
  if (type == 3){ // Farlie-Gumbel-Morgenstern copula
    result = pow(x, 1/paraX) * (1 - pow(x, 1/paraX)) * pow(y, 1/paraY) * (1 - pow(y, 1/paraY));
    res.push_back(result);
  }
  if (type == 4){ // Frank copula
    result = log(1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1)) / (a * a)
           - (-((exp(-a * pow(x, 1/paraX)) - 1) * pow(y, 1/paraY) * exp(-a * pow(y, 1/paraY)) / (exp(-a) - 1))
           -((exp(-a * pow(y, 1/paraY)) - 1) * pow(x, 1/paraX) * exp(-a * pow(x, 1/paraX)) / (exp(-a) - 1))
           +((exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1) * exp(-a) / pow((exp(-a) - 1), 2)))
           / (a * (1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1)));
    res.push_back(result);
  }
  if (type == 5) { // Gaussian copula
    /*Function func("graphi.norm");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A = wrap(a);
    NumericVector ret = func(X, PARAX, Y, PARAY, A); // Call function "graphi.norm" from R
    result = ret[0];
    res.push_back(result);
    * */
    result = grad(x, paraX, y, paraY, a); // Call "grad" function from Gradient.cpp
    res.push_back(result);
  }
  if (type == 6) { // Ali-Mikhail-Haq copula
    result = pow(x, 1/paraX) * pow(y, 1/paraY) * (1 - pow(y, 1/paraY)) * (1 - pow(x, 1/paraX))
           / pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 2);
    res.push_back(result);
  }
  if (type == 7) { // Joe copula
    result = - pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a)
           * (-log(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a)) / (a * a)
           + (-pow(1-pow(x, 1/paraX), a) * pow(1-pow(y, 1/paraY), a) * log((1-pow(x, 1/paraX)) * (1-pow(y, 1/paraY)))
           + pow(1 - pow(x, 1/paraX), a) * log(1 - pow(x, 1/paraX)) + pow(1 - pow(y, 1/paraY), a) * log(1 - pow(y, 1/paraY))) 
           / (a * (pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a))));
    res.push_back(result);
  }
  if (type == 8){ // Student copula
    //Function func1("graphi1.student");
    Function func2("graphi2.student");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A(2);
    A[0] = theta[0];
    A[1] = theta[1];
    //NumericVector ret1 = func1(X, PARAX, Y, PARAY, A); // Call function "graphi1.student" from R
    NumericVector ret2 = func2(X, PARAX, Y, PARAY, A); // Call function "graphi2.student" from R
    //double result1 = ret1[0];
    double result2 = ret2[0];
    //res.push_back(result1);
    res.push_back(result2);
  }
  return res;
}

// A gradient of partial derivative of one variable
vector<double>  graDxPhi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type){
  double result;
  vector<double> res;
  double a = theta[0];
  if (type == 0){ // Distribution given
    NumericVector ret = compute(x, paraX, y, paraY, a, f, g);
    Function dxgFunc("margin");
    NumericVector PARAX = wrap(paraX);
    NumericVector X = wrap(x);
    NumericVector deri = dxgFunc(PARAX, X, dxg);
    result = deri[0] * ret[0];
    res.push_back(result);
  }
  if (type == 1){ // Gumbel distribution
    result = - pow(x, -1/a - 1) * pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), 2 * a - 1)
         * (log(pow(paraX * x, -1/a) + pow(paraY * y, -1/a)) 
         + a * (log(paraX * x) * pow(paraX * x, -1/a) * pow(a,-2) + log(paraY * y) * pow(paraY * y, -1/a) * pow(a,-2) )
         / (pow(paraX * x, -1/a) + pow(paraY * y, -1/a)))
         * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX, -1/a)
         + pow(x, -1/a - 1) * pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a - 1)
         * (log(pow(paraX * x, -1/a) + pow(paraY * y, -1/a)) 
         + (a - 1) * (log(paraX * x) * pow(paraX * x, -1/a) * pow(a,-2) + log(paraY * y) * pow(paraY * y, -1/a) * pow(a,-2) )
         / (pow(paraX * x, -1/a) + pow(paraY * y, -1/a)))
         * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX, -1/a)
         + pow(x, -1/a - 1) * log(x) * pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a - 1) * pow(a,-2)
         * exp (- pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX, -1/a)
         + log(paraX) * pow(x, -1/a - 1) * pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a - 1) * pow(a,-2)
         * exp (- pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX, -1/a);
    res.push_back(result);
  }
  if (type == 2){ // Gumbel copula
    double a_tmp = a;
    a = (double) 1/a;
    result = pow(-log(x) / paraX, 1/a) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a)) * pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), 2 * a - 1)
           *   (log(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)) - ((pow(-log(y)/paraY, 1/a)*log(-log(y)/paraY)) + (pow(-log(x)/paraX, 1/a)*log(-log(x)/paraX)))
           /   (a * (pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)))) / (x * log(x))
           - pow(-log(x) / paraX, 1/a) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a)) * pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a - 1)
           *   (log(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)) - (a - 1) * ((pow(-log(y)/paraY, 1/a)*log(-log(y)/paraY)) + (pow(-log(x)/paraX, 1/a)*log(-log(x)/paraX)))
           /   (a * a * (pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)))) / (x * log(x))	
           +    pow(-log(x)/paraX, 1/a) * log(-log(x)/paraX) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a))
           *    pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a - 1) / (a * a * x * log(x));
    result = result * (-1/(a_tmp * a_tmp));
    res.push_back(result);
  }
  if (type == 3){ // Farlie-Gumbel-Morgenstern copula
    result = pow(x, 1/paraX - 1) * (1 - pow(x, 1/paraX)) * pow(y, 1/paraY) * (1 - pow(y, 1/paraY)) / paraX
           - pow(x, 2/paraX - 1) * pow(y, 1/paraY) * (1 - pow(y, 1/paraY)) / paraX;
    res.push_back(result);
  }
  if (type == 4){ // Frank copula
    result = - ((exp(-a * pow(y, 1/paraY)) - 1) * pow(x, 1/paraX - 1) * exp(-a * pow(x, 1/paraX)))
           * (-((exp(-a * pow(x, 1/paraX)) - 1) * pow(y, 1/paraY) * exp(-a * pow(y, 1/paraY)) / (exp(-a) - 1))
           -((exp(-a * pow(y, 1/paraY)) - 1) * pow(x, 1/paraX) * exp(-a * pow(x, 1/paraX)) / (exp(-a) - 1))
           +((exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1) * exp(-a) / pow((exp(-a) - 1), 2)))
           / (pow((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1)), 2) * paraX * (exp(-a) - 1))
           - (pow(x, 1/paraX - 1) * pow(y, 1/paraY) * exp(-a * pow(x, 1/paraX) - a * pow(y, 1/paraY)))
           / (((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1))) * paraX * (exp(-a) - 1))
           + (pow(x, 1/paraX - 1) * exp(-a * pow(x, 1/paraX) - a) * (exp(-a * pow(y, 1/paraY)) - 1))
           / (((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1))) * paraX * (exp(-a) - 1) * (exp(-a) - 1))
           - (pow(x, 2/paraX - 1) * exp(-a * pow(x, 1/paraX)) * (exp(-a * pow(y, 1/paraY)) - 1))
           / (((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1))) * paraX * (exp(-a) - 1));
    res.push_back(result);
  }
  if (type == 5) { // Gaussian copula
    /*Function func("gradxphi.norm");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A = wrap(a);
    NumericVector ret = func(X, PARAX, Y, PARAY, A); // Call function "gradxphi.norm" from R
    result = ret[0];
    res.push_back(result);
    * */
    double xx[n+1];
    double eps;
    eps = (double) pow(y,(double)1/paraY)/n;
    xx[0] = 0;
    for(int i = 1; i < n+1; i++) {
      xx[i] = xx[i-1] + eps;
    }
    // Compute the integral
    double summand = 0;
    for(int i=0; i<n; i++){
      double mean = (xx[i+1]+xx[i])/2;
      double val1 = Rf_qnorm5(mean, 0.0, 1.0, 1, 0);
      double val2 = Rf_qnorm5(pow(x,(double)1/paraX), 0.0, 1.0, 1, 0);
      double resultat = (-pow(a, 3) + val1 * val2 * pow(a, 2) + a * (1 - val1*val1 - val2*val2) + val1*val2)
                      * exp((pow(a * val1, 2) + pow(a * val2, 2) - 2 * a * val1 * val2)
                      / (-2 + 2 * a * a)) / pow((1 - a * a), 2.5);
      summand += resultat * eps;
    }
    summand *= pow(x, (double)1/paraX - 1) / paraX;
    res.push_back(summand);
  }
  if (type == 6) { // Ali-Mikhail-Haq copula
    result = - pow(x, 2/paraX - 1) * pow(y, 1/paraY) * (1 - pow(y, 1/paraY)) 
           / (paraX * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 2))
           + pow(x, 1/paraX - 1) * pow(y, 1/paraY) * (1 - pow(y, 1/paraY)) * (1 - pow(x, 1/paraX)) 
           / (paraX * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 2))
           - 2 * a * pow(x, 2/paraX - 1) * pow(y, 1/paraY) * pow((1 - pow(y, 1/paraY)), 2) * (1 - pow(x, 1/paraX)) 
           / (paraX * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 3));
    res.push_back(result);
  }
  if (type == 7) { // Joe copula
    result = - (pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * (pow(1 - pow(y, 1/paraY), a) - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 1)
           * (-log(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a)) / (a * a)
           + (1/a - 1) * (-pow(1-pow(x, 1/paraX), a) * pow(1-pow(y, 1/paraY), a) * log((1-pow(x, 1/paraX)) * (1-pow(y, 1/paraY)))
           + pow(1 - pow(x, 1/paraX), a) * log(1 - pow(x, 1/paraX)) + pow(1 - pow(y, 1/paraY), a) * log(1 - pow(y, 1/paraY))) 
           / (pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a)))) / paraX
           
           - pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 1)
           * (a * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * pow(1 - pow(y, 1/paraY), a) * 
           log((1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY)))
           + pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * pow(1 - pow(y, 1/paraY), a) 
           - a * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * log(1 - pow(x, 1/paraX))
           - pow(x, 1/paraX - 1) *  pow(1 - pow(x, 1/paraX), a - 1)) / (a * paraX)
           
           + pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * (pow(1 - pow(y, 1/paraY), a) - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 1) / (paraX * a);
    res.push_back(result);
  }
  if (type == 8){ // Student copula
    //Function func1("gradxphi1.student");
    Function func2("gradxphi2.student");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A(2);
    A[0] = theta[0];
    A[1] = theta[1];
    //NumericVector ret1 = func1(X, PARAX, Y, PARAY, A); // Call function "gradxphi1.student" from R
    NumericVector ret2 = func2(X, PARAX, Y, PARAY, A); // Call function "gradxphi2.student" from R
    //double result1 = ret1[0];
    double result2 = ret2[0];
    //res.push_back(result1);
    res.push_back(result2);
  }
  return res;
}

// Theta gradient of partial derivative of two variables
vector<double> graDxDyPhi(double x, double paraX, double y, double paraY, vector<double> theta, SEXP f, SEXP g, SEXP dxg, int type){
  double result;
  vector<double> res;
  double a = theta[0];
  if (type == 0){ // Distribution given
    NumericVector ret = compute(x, paraX, y, paraY, a, f, g);
    Function dxgFunc("margin");
    NumericVector PARAX = wrap(paraX);
    NumericVector X = wrap(x);
    NumericVector deri = dxgFunc(PARAX, X, dxg);
    NumericVector PARAY = wrap(paraY);
    NumericVector Y = wrap(y);
    NumericVector deri1 = dxgFunc(PARAY, Y, dxg);
    result = deri[0] * deri1[0] * ret[0];
    res.push_back(result);
  }
  if (type == 1){ // Gumbel distribution
    result = pow(x * y, -1/a -1) * pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), 2 * a - 2) 
         * (2 * log(pow(paraX * x, -1/a) + pow(paraY * y, -1/a)) + (2 * a - 2)
         * (log(paraX * x) * pow(a,-2) * pow(paraX * x, -1/a) + log(paraY * y) * pow(a,-2) * pow(paraY * y, -1/a))
         / (pow(paraX * x, -1/a) + pow(paraY * y, -1/a)))
         * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX * paraY, -1/a)
         - pow(x * y, -1/a -1) * pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), 3 * a - 2) 
         * (log(pow(paraX * x, -1/a) + pow(paraY * y, -1/a)) + a
         * (log(paraX * x) * pow(a,-2) * pow(paraX * x, -1/a) + log(paraY * y) * pow(a,-2) * pow(paraY * y, -1/a))
         / (pow(paraX * x, -1/a) + pow(paraY * y, -1/a)))
         * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX * paraY, -1/a)
         + (1 - 1/a) * pow(x * y, -1/a -1) * pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), 2 * a - 2) 
         * (log(pow(paraX * x, -1/a) + pow(paraY * y, -1/a)) + a
         * (log(paraX * x) * pow(a,-2) * pow(paraX * x, -1/a) + log(paraY * y) * pow(a,-2) * pow(paraY * y, -1/a))
         / (pow(paraX * x, -1/a) + pow(paraY * y, -1/a)))
         * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX * paraY, -1/a)
         - (1 - 1/a) * pow(x * y, -1/a -1) * pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a - 2) 
         * (log(pow(paraX * x, -1/a) + pow(paraY * y, -1/a)) + (a - 2)
         * (log(paraX * x) * pow(a,-2) * pow(paraX * x, -1/a) + log(paraY * y) * pow(a,-2) * pow(paraY * y, -1/a))
         / (pow(paraX * x, -1/a) + pow(paraY * y, -1/a)))
         * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX * paraY, -1/a)
         + pow(x * y, -1/a - 1) * pow (pow(paraX * x, -1/a) + pow(paraY * y, -1/a), 2 * a - 2) * log(x * y)
         * pow(a, -2) * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX * paraY, -1/a)
         - (1 - 1/a) * pow(x * y, -1/a - 1) * pow (pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a - 2) * log(x * y)
         * pow(a, -2) * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX * paraY, -1/a)
         + log(paraX * paraY) * pow(x * y, -1/a - 1) * pow (pow(paraX * x, -1/a) + pow(paraY * y, -1/a), 2 * a - 2)
         * pow(a, -2) * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX * paraY, -1/a)
         - (1 - 1/a) * pow(x * y, -1/a - 1) * pow (pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a - 2) * log(paraX * paraY)
         * pow(a, -2) * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX * paraY, -1/a)
         - pow(x * y, -1/a - 1) * pow (pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a - 2)
         * pow(a, -2) * exp(-pow(pow(paraX * x, -1/a) + pow(paraY * y, -1/a), a)) * pow(paraX * paraY, -1/a);
    res.push_back(result);
  }
  if (type == 2){ // Gumbel copula
    double a_tmp = a;
    a = (double) 1/a;
    result = pow(-log(x) / paraX, 1/a) * pow(-log(y) / paraY, 1/a) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a)) * pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), 2 * a - 2)
           *   (2 * log(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)) - (2 * a - 2)* ((pow(-log(y)/paraY, 1/a)*log(-log(y)/paraY)) + (pow(-log(x)/paraX, 1/a)*log(-log(x)/paraX)))
           /   (a * a * (pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)))) / (x * log(x) * y * log(y))
           - pow(-log(x) / paraX, 1/a) * pow(-log(y) / paraY, 1/a) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a)) * pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), 3 * a - 2)
           *   (log(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)) - (a)* ((pow(-log(y)/paraY, 1/a)*log(-log(y)/paraY)) + (pow(-log(x)/paraX, 1/a)*log(-log(x)/paraX)))
           /   (a * a * (pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)))) / (x * log(x) * y * log(y))
           + (a - 1) * pow(-log(x) / paraX, 1/a) * pow(-log(y) / paraY, 1/a) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a)) * pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), 2 * a - 2)
           *   (log(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)) - (a)* ((pow(-log(y)/paraY, 1/a)*log(-log(y)/paraY)) + (pow(-log(x)/paraX, 1/a)*log(-log(x)/paraX)))
           /   (a * a * (pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)))) / (a * x * log(x) * y * log(y))
           - (a - 1) * pow(-log(x) / paraX, 1/a) * pow(-log(y) / paraY, 1/a) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a)) * pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a - 2)
           *   (log(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)) - (a - 2)* ((pow(-log(y)/paraY, 1/a)*log(-log(y)/paraY)) + (pow(-log(x)/paraX, 1/a)*log(-log(x)/paraX)))
           /   (a * a * (pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a)))) / (a * x * log(x) * y * log(y))
           -    pow(-log(x)/paraX, 1/a) * pow(-log(y)/paraY, 1/a) * (log(-log(y)/paraY)+log(-log(x)/paraX)) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a))
           *    pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), 2 * a - 2) / (a * a * x * log(x) * y * log(y))
           +    (a - 1) * pow(-log(x)/paraX, 1/a) * pow(-log(y)/paraY, 1/a) * (log(-log(y)/paraY)+log(-log(x)/paraX)) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a))
           *    pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a - 2) / (a * a * a * x * log(x) * y * log(y))
           -    pow(-log(x)/paraX, 1/a) * pow(-log(y)/paraY, 1/a) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a))
           *    pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a - 2) / (a * x * log(x) * y * log(y))
           +    (a - 1) * pow(-log(x)/paraX, 1/a) * pow(-log(y)/paraY, 1/a) * exp(-pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a))
           *    pow(pow(-log(x) / paraX, 1/a) + pow(-log(y) / paraY, 1/a), a - 2) / (a * a * x * log(x) * y * log(y));
    result = result * (-1/(a_tmp * a_tmp));
    res.push_back(result);
  }
  if (type == 3){ // Farlie-Gumbel-Morgenstern copula
    result = pow(x, 2/paraX - 1) * pow(y, 2/paraY - 1) / (paraX * paraY)
           - pow(x, 1/paraX - 1) * (1 - pow(x, 1/paraX)) * pow(y, 2/paraY - 1) / (paraX * paraY)
           - pow(y, 1/paraY - 1) * (1 - pow(y, 1/paraY)) * pow(x, 2/paraX - 1) / (paraX * paraY)
           + pow(x, 1/paraX - 1) * pow(y, 1/paraY -  1) * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY)) / (paraX * paraY);
    res.push_back(result);
  }
  if (type == 4){ // Frank copula
    result = (a * pow(x, 1/paraX - 1) * pow(y, 1/paraY - 1) * exp(-a * pow(x, 1/paraX) - a * pow(y, 1/paraY) - a) * (exp(-a * pow(y, 1/paraY)) - 1) * (exp(-a * pow(x, 1/paraX)) - 1))
           / (pow(((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1))), 2) * paraX * paraY * pow((exp(-a) - 1), 3))
           - (a * pow(x, 2/paraX - 1) * pow(y, 1/paraY - 1) * exp(-a * pow(x, 1/paraX) - a * pow(y, 1/paraY)) * (exp(-a * pow(y, 1/paraY)) - 1) * (exp(-a * pow(x, 1/paraX)) - 1))
           / (pow(((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1))), 2) * paraX * paraY * pow((exp(-a) - 1), 2))
           - (pow(x, 1/paraX - 1) * pow(y, 1/paraY - 1) * exp(-a * pow(x, 1/paraX) - a * pow(y, 1/paraY)))
           / (((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1))) * paraX * paraY * (exp(-a) - 1))
           - (a * pow(x, 1/paraX - 1) * pow(y, 2/paraY - 1) * exp(-a * pow(x, 1/paraX) - 2 * a * pow(y, 1/paraY)) * (exp(-a * pow(x, 1/paraX)) - 1))
           / (pow(((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1))), 2) * paraX * paraY * pow((exp(-a) - 1), 2))           
           + (a * pow(x, 1/paraX - 1) * pow(y, 2/paraY - 1) * exp(-a * pow(x, 1/paraX) - a * pow(y, 1/paraY)))
           / (((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1))) * paraX * paraY * (exp(-a) - 1))
           - (a * pow(x, 1/paraX - 1) * pow(y, 1/paraY - 1) * exp(-a * pow(x, 1/paraX) - a * pow(y, 1/paraY) - a))
           / (((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1))) * paraX * paraY * pow((exp(-a) - 1), 2))
           + (a * pow(x, 2/paraX - 1) * pow(y, 1/paraY - 1) * exp(-a * pow(x, 1/paraX) - a * pow(y, 1/paraY)))
           / (((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1))) * paraX * paraY * (exp(-a) - 1))           
           - ((exp(-a * pow(y, 1/paraY)) - 1) * pow(x, 1/paraX - 1) * exp(-a * pow(x, 1/paraX)))
           * ((a * (exp(-a * pow(x, 1/paraX)) - 1) * pow(y, 2/paraY - 1) * exp(-a * pow(y, 1/paraY)) / (paraY * (exp(-a) - 1)))
           + (a * pow(y, 1/paraY - 1) * pow(x, 1/paraX) * exp(-a * pow(y, 1/paraY) - a * pow(x, 1/paraX)) / (paraY * (exp(-a) - 1)))
           - (a * (exp(-a * pow(x, 1/paraX)) - 1) * pow(y, 1/paraY - 1) * exp(-a * pow(y, 1/paraY) - a) / (paraY * pow((exp(-a) - 1), 2)))
           - ((exp(-a * pow(x, 1/paraX)) - 1) * pow(y, 1/paraY - 1) * exp(-a * pow(y, 1/paraY)) / (paraY * (exp(-a) - 1))))
           / (pow((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1)), 2) * paraX * (exp(-a) - 1))           
           - (2 * a * (exp(-a * pow(y, 1/paraY)) - 1) * (exp(-a * pow(x, 1/paraX)) - 1) * pow(x, 1/paraX - 1) * pow(y, 1/paraY - 1) * exp(-a * pow(x, 1/paraX) - a * pow(y, 1/paraY)))
           * (-((exp(-a * pow(x, 1/paraX)) - 1) * pow(y, 1/paraY) * exp(-a * pow(y, 1/paraY)) / (exp(-a) - 1))
           -((exp(-a * pow(y, 1/paraY)) - 1) * pow(x, 1/paraX) * exp(-a * pow(x, 1/paraX)) / (exp(-a) - 1))
           +((exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1) * exp(-a) / pow((exp(-a) - 1), 2)))
           / (pow((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1)), 3) * paraX * paraY * pow((exp(-a) - 1), 2))           
           + (a * pow(x, 1/paraX - 1) * pow(y, 1/paraY - 1) * exp(-a * pow(x, 1/paraX) - a * pow(y, 1/paraY)))
           * (-((exp(-a * pow(x, 1/paraX)) - 1) * pow(y, 1/paraY) * exp(-a * pow(y, 1/paraY)) / (exp(-a) - 1))
           -((exp(-a * pow(y, 1/paraY)) - 1) * pow(x, 1/paraX) * exp(-a * pow(x, 1/paraX)) / (exp(-a) - 1))
           +((exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1) * exp(-a) / pow((exp(-a) - 1), 2)))
           / (pow((1 + (exp(-a * pow(x, 1/paraX)) - 1) * (exp(-a * pow(y, 1/paraY)) - 1)/(exp(-a) - 1)), 2) * paraX * paraY * (exp(-a) - 1));
    res.push_back(result);
  }
  if (type == 5) { // Gaussian copula
    /*Function func("gradxdyphi.norm");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A = wrap(a);
    NumericVector ret = func(X, PARAX, Y, PARAY, A); // Call function "gradxdyphi.norm" from R
    result = ret[0];
    res.push_back(result);
    * */    
    double val1 = Rf_qnorm5(pow(x, 1/paraX), 0.0, 1.0, 1, 0);
    double val2 = Rf_qnorm5(pow(y, 1/paraY), 0.0, 1.0, 1, 0);
    result = pow(x, 1/paraX - 1) * pow (y, 1/paraY - 1) / (paraX * paraY);
    result *= (-pow(a, 3) + val1 * val2 * pow(a, 2) + a * (1 - val1*val1 - val2*val2) + val1*val2)
           * exp((pow(a * val1, 2) + pow(a * val2, 2) - 2 * a * val1 * val2)
           / (-2 + 2 * a * a)) / pow((1 - a * a), 2.5);
    res.push_back(result);
  }
  if (type == 6) { // Ali-Mikhail-Haq copula
    result = - pow(x, 2/paraX - 1) * pow(y, 1/paraY - 1) * (1 - pow(y, 1/paraY)) 
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 2))
           + pow(x, 1/paraX - 1) * pow(y, 1/paraY - 1) * (1 - pow(y, 1/paraY)) * (1 - pow(x, 1/paraX)) 
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 2))
           + 6 * a * pow(x, 2/paraX - 1) * pow(y, 2/paraY - 1) * (1 - pow(y, 1/paraY)) * (1 - pow(x, 1/paraX)) 
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 3))
           - 2 * a * pow(x, 1/paraX - 1) * pow(y, 2/paraY - 1) * (1 - pow(y, 1/paraY)) * pow((1 - pow(x, 1/paraX)), 2)
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 3))
           - 2 * a * pow(x, 2/paraX - 1) * pow(y, 1/paraY - 1) * (1 - pow(x, 1/paraX)) * pow((1 - pow(y, 1/paraY)), 2)
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 3))
           + 6 * a * a * pow(x, 2/paraX - 1) * pow(y, 2/paraY - 1) * pow((1 - pow(y, 1/paraY)), 2) * pow((1 - pow(x, 1/paraX)), 2) 
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))),4))
           + pow(x, 2/paraX - 1) * pow(y, 2/paraY - 1)
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 2))
           - pow(x, 1/paraX - 1) * pow(y, 2/paraY - 1) * (1 - pow(x, 1/paraX)) 
           / (paraX * paraY * pow((1 - a * (1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY))), 2));
    res.push_back(result);
  }
  if (type == 7) { // Joe copula
    result = a * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 1)
           * ((1/a - 1) * (-pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a) * log((1 - pow(y, 1/paraY)) * (1 - pow(x, 1/paraX)))
           + pow(1 - pow(y, 1/paraY), a) * log(1 - pow(y, 1/paraY)) + pow(1 - pow(x, 1/paraX), a) * log(1 - pow(x, 1/paraX))) 
           / (pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a))
           - log(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a)) / (a * a)
           ) / (paraX * paraY)         
           
           - (1 - a) * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1)
           * (pow(1 - pow(x, 1/paraX), a) - 1) * (pow(1 - pow(y, 1/paraY), a) - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 2)
           * ((1/a - 2) * (-pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a) * log((1 - pow(y, 1/paraY)) * (1 - pow(x, 1/paraX)))
           + pow(1 - pow(y, 1/paraY), a) * log(1 - pow(y, 1/paraY)) + pow(1 - pow(x, 1/paraX), a) * log(1 - pow(x, 1/paraX))) 
           / (pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a))
           - log(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a)
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a)) / (a * a)
           ) / (paraX * paraY)
          
           - (1 - a) * pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1) * (pow(1 - pow(x, 1/paraX), a) - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 2)
           * (a * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * pow(1 - pow(y, 1/paraY), a) * 
           log((1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY)))
           + pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * pow(1 - pow(y, 1/paraY), a) 
           - a * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * log(1 - pow(x, 1/paraX))
           - pow(x, 1/paraX - 1) *  pow(1 - pow(x, 1/paraX), a - 1)) / (a * paraX * paraY)
           
           - (1 - a) * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * (pow(1 - pow(y, 1/paraY), a) - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 2)
           * (a * pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1) * pow(1 - pow(x, 1/paraX), a) * 
           log((1 - pow(x, 1/paraX)) * (1 - pow(y, 1/paraY)))
           + pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1) * pow(1 - pow(x, 1/paraX), a) 
           - a * pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1) * log(1 - pow(y, 1/paraY))
           - pow(y, 1/paraY - 1) *  pow(1 - pow(y, 1/paraY), a - 1)) / (a * paraX * paraY)

           + a * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 1) 
           * log((1 - pow(y, 1/paraY)) * (1 - pow(x, 1/paraX))) / (paraX * paraY)           
           
           + pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1) * pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 1) / (paraX * paraY)
           
           + (1 - a) * pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1) * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1)
           * (pow(1 - pow(x, 1/paraX), a) - 1) * (pow(1 - pow(y, 1/paraY), a) - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 2) / (paraX * paraY * a)
           
           + pow(y, 1/paraY - 1) * pow(1 - pow(y, 1/paraY), a - 1) * pow(x, 1/paraX - 1) * pow(1 - pow(x, 1/paraX), a - 1)
           * (pow(1 - pow(x, 1/paraX), a) - 1) * (pow(1 - pow(y, 1/paraY), a) - 1)
           * pow(pow(1 - pow(x, 1/paraX), a) + pow(1 - pow(y, 1/paraY), a) 
           - pow(1 - pow(x, 1/paraX), a) * pow(1 - pow(y, 1/paraY), a), 1/a - 2) / (paraX * paraY * a);
    res.push_back(result);
  }
  if (type == 8){ // Student copula
    //Function func1("gradxdyphi1.student");
    Function func2("gradxdyphi2.student");
    NumericVector PARAX = wrap(paraX);
    NumericVector PARAY = wrap(paraY);
    NumericVector X = wrap(x); // Convert from double in C++ to SEXP in Rcpp
    NumericVector Y = wrap(y);
    NumericVector A(2);
    A[0] = theta[0];
    A[1] = theta[1];
    //NumericVector ret1 = func1(X, PARAX, Y, PARAY, A); // Call function "gradxdyphi1.student" from R
    NumericVector ret2 = func2(X, PARAX, Y, PARAY, A); // Call function "gradxdyphi2.student" from R
    //double result1 = ret1[0];
    double result2 = ret2[0];
    //res.push_back(result1);
    res.push_back(result2);
  }
  return res;
}
