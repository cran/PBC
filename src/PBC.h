/*
  Copyright (C) 2013 Pham Van Trung and Mazo Gildas

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 3 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Derivative-sum-product (DSP) algorithm (by Huang and Frey)
 * and gradient-derivative-product (GDP) (by Huang and Nebojsa Jojic)
 * 
 * References: 
 * 1. Cumulative distribution networks: Inference, estimation and applications of graphical models 
 * for cumulative distribution functions. Jim C. Huang (2009). Ph.D. Thesis. 
 * 2. Maximum-likelihood learning of cumulative distribution functions on graphs. 
 * Jim C. Huang and Nebojsa Jojic (2010). Microsoft Research. Redmond, WA 98052.
 */
 
#ifndef PBC_H
#define PBC_H

#include <vector>

#include <Rcpp.h>
#include <math.h>

using namespace std;
using namespace Rcpp;

class PBC {    
  public :
    void setOut(SEXP m_out);
    void setDxg(SEXP m_dxg);
    void setMargins(SEXP m_g);
    void setType(SEXP m_type);
    void setNIteration(SEXP m_nIteration);
    void setFunction(SEXP m_f);
    void setDxf(SEXP m_dxf);
	void setDxdyf(SEXP m_dxdyf);
	void setGraf(SEXP m_graf);
	void setGradxf(SEXP m_gradxf);
	void setGradxdyf(SEXP m_gradxdyf);
    void setRoot(SEXP m_root);
    void setVector(SEXP m_x);
    void setBinMat(SEXP m_binMat, SEXP m_x);
    void setTheta(SEXP m_theta);
    vector <double> multiple(vector<double> vect, double x);
    vector <double> add(vector<double> vect1, vector<double> vect2);
    vector <double> addNumberVector(double x, int pos, vector<double> vect);
    void algoInit(vector <bool> isLeaf_v);
    void parametersInit(SEXP m_binMat, SEXP m_theta, SEXP m_x, SEXP m_root, SEXP m_f, SEXP m_nIteration, 
                        SEXP m_dxf, SEXP m_dxdyf, SEXP m_graf, SEXP m_gradxf, SEXP m_gradxdyf, SEXP m_type, SEXP m_g, SEXP m_dxg, SEXP out);
    NumericVector algo();
    PBC () {}
    
  private:
    SEXP g_; // Margin function
    SEXP dxg_; // Derivative of margin function
    SEXP f_; // Distribution function
    SEXP dxf_; // Partial derivative of one variable
    SEXP dxdyf_; // Partial derivative of two variables
    SEXP graf_; // Gradient of theta
    SEXP gradxf_; // Theta gradient of partial derivative of one variable
    SEXP gradxdyf_; // Theta gradient of partial derivative of two variables
    int root_; // Root number
    int nIteration_; // Number of interation
    vector<double> vect_; // Data
    vector< vector<int> > binMat_; // Relation between functions and variables
    vector< vector<double> > fv_mu_; // From functions to variables - mu
    vector< vector<double> > fv_la_; // From functions to variables - lamda
    vector< vector<double> > vf_mu_; // From variables to functions - mu
    vector< vector<double> > vf_la_; // From variables to functions - lamda
    vector< vector< vector<double> > > gra_fv_mu_; // Gradient from functions to variables - mu
    vector< vector< vector<double> > > gra_fv_la_; // Gradient from functions to variables - lamda
    vector< vector< vector<double> > > gra_vf_mu_; // Gradient from variables to functions - mu
    vector< vector< vector<double> > > gra_vf_la_; // Gradient from variables to functions - lamda 
    vector<double> theta_; // Function parameter theta
    int type_; // Type of distribution
    int out_; // Display resultats

};

#endif // PBC_H
