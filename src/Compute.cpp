/*
 * Compute.cpp : Tranfer resultats from C++ to R 
 */

#include <Rcpp.h>

#include "PBC.h"
#include "Derivatives.h"

using namespace Rcpp;

RcppExport SEXP pbc(SEXP m_binMat, SEXP m_theta, SEXP m_x, SEXP m_root, SEXP m_f, SEXP m_nIteration, SEXP m_dxf, SEXP m_dxdyf, 
                    SEXP m_graf, SEXP m_gradxf, SEXP m_gradxdyf, SEXP m_type, SEXP m_g, SEXP m_dxg, SEXP out) {
  PBC pbc;
  pbc.parametersInit(m_binMat, m_theta, m_x, m_root, m_f, m_nIteration, m_dxf, m_dxdyf, m_graf, m_gradxf, m_gradxdyf, m_type, m_g, m_dxg, out); // Parameters transmission 
  NumericVector resultat;
  resultat = pbc.algo(); // Call message-passing algorithm	  
  return(resultat);
}
