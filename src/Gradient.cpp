/*
 * Richardson Method
 * Applying Richardson Extrapolation to improve the accuracy of 
 * the first and second order derivatives. The algorithm as follows:
 * --  For each column of the derivative matrix a,
 * say, A1, A2, ..., Ar, by Richardson Extrapolation, to calculate a
 * new sequence of approximations B1, B2, ..., Br used the formula
 * 
 * B(i) =( A(i+1)*4^m - A(i) ) / (4^m - 1) ,  i=1,2,...,r-m
 * 
 * N.B. This formula assumes v=2.
 * -- Initially m is taken as 1  and then the process is repeated 
 * restarting with the latest improved values and increasing the 
 * value of m by one each until m equals r-1
 * 
*/

/*
 * 
 * Modified and translated R code from 'numDeriv' package (by Paul Gilbert) to C++ code (by Pham Van Trung)
 * 
 */

#include "Gradient.h"
#include "miwa.h"

using namespace std;

// Compute numerically gradient of function
double grad(double u1, double n1, double u2, double n2, double x){
  double result;
  double zeroTol = 178102.9;
  double eps = 0.0001;
  double d = 0.0001;
  int r = 4;
  int v = 2;
  double* a = new double[r];
  double h = 0.0;
  if (fabs(x) < zeroTol)
    h = fabs(d * x) + eps;
  else
    h = fabs(d * x);
  for (int i = 0; i < r; i++) {
    a[i] = (calcul(u1, n1, u2, n2, x + h) - calcul(u1, n1, u2, n2, x - h)) / (2 * h);
    h = (double) h / v;
  }
  double* tmp = a;  
  for (int i = 0; i < (r-1); i++) {
    for (int j = 0; j < r - i - 1; j++)
      tmp[j] = (a[j + 1] * pow( 4.0, (double) i + 1) - a[j]) / (pow(4.0, (double) i + 1) - 1);
    a = tmp;  
  }
  result = a[0];
  return result;
}

