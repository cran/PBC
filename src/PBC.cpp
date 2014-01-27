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
 
#include <iostream>
#include <algorithm>

#include "Derivatives.h"
#include "PBC.h"

// Setter functions

void PBC::setOut(SEXP m_out) {
  out_ = as<int>(m_out);
}



void PBC::setDxg(SEXP m_dxg) {
  dxg_ = m_dxg;
}



void PBC::setMargins(SEXP m_g) {
  g_ = m_g;
}



void PBC::setType(SEXP m_type) {
  type_ = as<int>(m_type);
}



void PBC::setNIteration(SEXP m_nIteration) {
  nIteration_ = as<int>(m_nIteration);
}



void PBC::setFunction(SEXP m_f) {
  f_ = m_f;
}



void PBC::setDxf(SEXP m_dxf) {
  dxf_ = m_dxf;
}



void PBC::setDxdyf(SEXP m_dxdyf) {
  dxdyf_ = m_dxdyf;
}



void PBC::setGraf(SEXP m_graf) {
  graf_ = m_graf;
}

void PBC::setGradxf(SEXP m_gradxf) {
  gradxf_ = m_gradxf;
}



void PBC::setGradxdyf(SEXP m_gradxdyf) {
  gradxdyf_ = m_gradxdyf;
}

void PBC::setRoot(SEXP m_root) {
  root_ = as<int>(m_root) - 1;
}



void PBC::setVector(SEXP m_x) {
  NumericVector x(m_x);
  vect_.resize(x.size());
  for (int i = 0; i < x.size(); i++)
    vect_[i] = x[i];
}



void PBC::setBinMat(SEXP m_binMat, SEXP m_x) {
  NumericVector x(m_x);
  NumericVector _binMat(m_binMat);
  binMat_.resize(x.size());
  for (unsigned i = 0; i < binMat_.size(); i++)
    binMat_[i].resize(x.size() - 1);
  for (int i = 0; i < _binMat.size(); i++)
    binMat_[i % x.size()][i / x.size()] = _binMat[i];
}



void PBC::setTheta(SEXP m_theta) {
  NumericVector a(m_theta);
  theta_.resize(a.size());
  for (int i = 0; i < a.size(); i++)
    theta_[i] = a[i];  
}



// Multiple of number and vector
vector <double> PBC::multiple(vector<double> vect, double x) {
  for (unsigned i = 0; i < vect.size(); i++)
    vect[i] *= x;
  return vect;
}



// Sum of two vectors
vector <double> PBC::add(vector<double> vect1, vector<double> vect2) {
  for (unsigned i = 0; i < vect1.size(); i++)
    vect1[i] += vect2[i];
  return vect1;
}



// Add a number to an element of vector
vector <double> PBC::addNumberVector(double x, int pos, vector<double> vect) {
  vect[pos] += x;
  return vect;
}



// Initialization of algorithm
void PBC::algoInit(vector <bool> isLeaf_v) {
  // Initialisation for mu and lamda functions
  vf_mu_.resize(binMat_.size());
  gra_vf_mu_.resize(binMat_.size());
  vector <double> tmp(binMat_[0].size());
  vector < vector <double> > tmp1(binMat_[0].size());
  vector <double> tmp2(theta_.size());
  for (unsigned i = 0; i < binMat_.size(); i++) {
    vf_mu_[i] = tmp;
    for (unsigned j = 0; j < binMat_[i].size(); j++) {	
	  tmp1[j] = tmp2;
	}
	gra_vf_mu_[i] = tmp1;
  }
  vf_la_ = vf_mu_;
  fv_mu_ = vf_mu_;
  fv_la_ = vf_mu_;
  gra_vf_la_ = gra_vf_mu_;
  gra_fv_mu_ = gra_vf_mu_;
  gra_fv_la_ = gra_vf_mu_;
  vector<double>().swap(tmp); // Free memory
  vector < vector <double> >().swap(tmp1); // Free memory
  vector<double>().swap(tmp2); // Free memory
  for (unsigned i = 0; i < binMat_.size(); i++)
    for (unsigned j = 0; j < binMat_[i].size(); j++)
      if ((isLeaf_v[i])&&(binMat_[i][j] == 1))
        vf_mu_[i][j] = 1.0;
}



// Initialization of parameters from R
void PBC::parametersInit(SEXP m_binMat, SEXP m_theta, SEXP m_x, SEXP m_root, SEXP m_f, SEXP m_nIteration, 
                         SEXP m_dxf, SEXP m_dxdyf, SEXP m_graf, SEXP m_gradxf, SEXP m_gradxdyf, SEXP m_type, SEXP m_g, SEXP m_dxg, SEXP m_out) {
  setType(m_type);
  setNIteration(m_nIteration);
  setFunction(m_f);
  setDxf(m_dxf);
  setDxdyf(m_dxdyf);
  setGraf(m_graf);
  setGradxf(m_gradxf);
  setGradxdyf(m_gradxdyf);
  setVector(m_x);  
  setTheta(m_theta);  
  setRoot(m_root);  
  setBinMat(m_binMat, m_x);
  setMargins(m_g);
  setDxg(m_dxg);
  setOut(m_out);
}



NumericVector PBC::algo() {
  unsigned int nFunc = binMat_[0].size(); // Number of function
  unsigned int nVar = binMat_.size(); // Number of variable (nFunc = nVar + 1)
  unsigned int nPara = theta_.size(); // Total number of parameters
  vector<int> nbNeighbor(nVar); // Number neighbor of variable
  unsigned nbTheta = (int)(nPara/nFunc); // Number of parameters for each function
  
  vector< vector <double> > a(nFunc); // Vector of parameters
  for (int unsigned i = 0; i < nFunc; i++) {
    vector<double> tmp(nbTheta);
    for (int unsigned j = 0; j < nbTheta; j++)
      tmp[j] = theta_[nbTheta * i + j];
	a[i] = tmp;
	vector<double>().swap(tmp);
  }  
  vector <bool> isLeaf_v(nVar);
  
  // Find the number neighbor of variable and leaf variables
  for (unsigned i = 0; i < nVar; i++) {
    int tmp = 0;
    for (unsigned j = 0; j < nFunc; j++)
      tmp += binMat_[i][j];
    nbNeighbor[i] = tmp;
    if (tmp > 1)
      isLeaf_v[i] = false;
    else
      isLeaf_v[i] = true;
  }

  algoInit(isLeaf_v); // Initialisation of algorithm
  
  NumericVector ret(nPara+1); // Results (P and gradient of P)
  bool isStop = false; // For stopping algorithm
  while (!isStop){  
    // Messages from functions to variables
    for (int unsigned i = 0; i < nFunc; i++) {
	  vector <int> neiVar; // Neighbor variables of function
	  for (int unsigned j = 0; j < nVar; j++)
	    if (binMat_[j][i] == 1)
	      neiVar.push_back(j);
	    
	  // Update mu and lamda functions from function to variable  
	  fv_mu_[neiVar[1]][i] = dxPhi(vect_[neiVar[0]], nbNeighbor[neiVar[0]], vect_[neiVar[1]], nbNeighbor[neiVar[1]], a[i], dxf_, g_, dxg_, type_) * vf_mu_[neiVar[0]][i]
	              + Phi(vect_[neiVar[0]], nbNeighbor[neiVar[0]], vect_[neiVar[1]], nbNeighbor[neiVar[1]], a[i], f_, g_, dxg_, type_) * vf_la_[neiVar[0]][i];
	  fv_mu_[neiVar[0]][i] = dxPhi(vect_[neiVar[1]], nbNeighbor[neiVar[1]], vect_[neiVar[0]], nbNeighbor[neiVar[0]], a[i], dxf_, g_, dxg_, type_) * vf_mu_[neiVar[1]][i]
	              + Phi(vect_[neiVar[1]], nbNeighbor[neiVar[1]], vect_[neiVar[0]], nbNeighbor[neiVar[0]], a[i], f_, g_, dxg_, type_) * vf_la_[neiVar[1]][i];               
	  fv_la_[neiVar[1]][i] = dxdyPhi(vect_[neiVar[0]], nbNeighbor[neiVar[0]], vect_[neiVar[1]], nbNeighbor[neiVar[1]], a[i], dxdyf_, g_, dxg_, type_) * vf_mu_[neiVar[0]][i]
	              + dxPhi(vect_[neiVar[1]], nbNeighbor[neiVar[1]], vect_[neiVar[0]], nbNeighbor[neiVar[0]], a[i], dxf_, g_, dxg_, type_) * vf_la_[neiVar[0]][i];
	  fv_la_[neiVar[0]][i] = dxdyPhi(vect_[neiVar[1]], nbNeighbor[neiVar[1]], vect_[neiVar[0]], nbNeighbor[neiVar[0]], a[i], dxdyf_, g_, dxg_, type_) * vf_mu_[neiVar[1]][i]
	              + dxPhi(vect_[neiVar[0]], nbNeighbor[neiVar[0]], vect_[neiVar[1]], nbNeighbor[neiVar[1]], a[i], dxf_, g_, dxg_, type_) * vf_la_[neiVar[1]][i];
	  
	  if (out_ > 0){
	    vector <double> grPhi(nPara), grDxPhi1(nPara), grDxPhi2(nPara), grDxDyPhi(nPara);
	    for (unsigned int ii = 0; ii < nbTheta; ii++){
	      grPhi = addNumberVector(graPhi(vect_[neiVar[0]], nbNeighbor[neiVar[0]], vect_[neiVar[1]], nbNeighbor[neiVar[1]], a[i], graf_, g_, dxg_, type_)[ii], i * nbTheta + ii, grPhi);
	      grDxPhi1 = addNumberVector(graDxPhi(vect_[neiVar[0]], nbNeighbor[neiVar[0]], vect_[neiVar[1]], nbNeighbor[neiVar[1]], a[i], gradxf_, g_, dxg_, type_)[ii], i * nbTheta + ii, grDxPhi1);
	      grDxPhi2 = addNumberVector(graDxPhi(vect_[neiVar[1]], nbNeighbor[neiVar[1]], vect_[neiVar[0]], nbNeighbor[neiVar[0]], a[i], gradxf_, g_, dxg_, type_)[ii], i * nbTheta + ii, grDxPhi2);
	      grDxDyPhi = addNumberVector(graDxDyPhi(vect_[neiVar[1]], nbNeighbor[neiVar[1]], vect_[neiVar[0]], nbNeighbor[neiVar[0]], a[i], gradxdyf_, g_, dxg_, type_)[ii], i * nbTheta + ii, grDxDyPhi);
	    }
	    // Update gradient
	    gra_fv_mu_[neiVar[1]][i] = add(add(add((multiple(grDxPhi1,vf_mu_[neiVar[0]][i])),
	                              (multiple(gra_vf_mu_[neiVar[0]][i],dxPhi(vect_[neiVar[0]], nbNeighbor[neiVar[0]], vect_[neiVar[1]], nbNeighbor[neiVar[1]], a[i], dxf_, g_, dxg_, type_)))),
	                              (multiple(grPhi,vf_la_[neiVar[0]][i]))),
			   					  (multiple(gra_vf_la_[neiVar[0]][i],Phi(vect_[neiVar[0]], nbNeighbor[neiVar[0]], vect_[neiVar[1]], nbNeighbor[neiVar[1]], a[i], f_, g_, dxg_, type_))));
	    gra_fv_mu_[neiVar[0]][i] = add(add(add((multiple(grDxPhi2,vf_mu_[neiVar[1]][i])),
	                              (multiple(gra_vf_mu_[neiVar[1]][i],dxPhi(vect_[neiVar[1]], nbNeighbor[neiVar[1]], vect_[neiVar[0]], nbNeighbor[neiVar[0]], a[i], dxf_, g_, dxg_, type_)))),
	                              (multiple(grPhi,vf_la_[neiVar[1]][i]))),
	                              (multiple(gra_vf_la_[neiVar[1]][i],Phi(vect_[neiVar[1]], nbNeighbor[neiVar[1]], vect_[neiVar[0]], nbNeighbor[neiVar[0]], a[i], f_, g_, dxg_, type_))));
	                          
	    gra_fv_la_[neiVar[1]][i] = add(add(add((multiple(grDxDyPhi,vf_mu_[neiVar[0]][i])),
	                              (multiple(gra_vf_mu_[neiVar[0]][i],dxdyPhi(vect_[neiVar[0]], nbNeighbor[neiVar[0]], vect_[neiVar[1]], nbNeighbor[neiVar[1]], a[i], dxdyf_, g_, dxg_, type_)))),
	                              (multiple(grDxPhi2,vf_la_[neiVar[0]][i]))),
	                              (multiple(gra_vf_la_[neiVar[0]][i],dxPhi(vect_[neiVar[1]], nbNeighbor[neiVar[1]], vect_[neiVar[0]], nbNeighbor[neiVar[0]], a[i], dxf_, g_, dxg_, type_))));
	                          
	    gra_fv_la_[neiVar[0]][i] = add(add(add((multiple(grDxDyPhi,vf_mu_[neiVar[1]][i])),
	                              (multiple(gra_vf_mu_[neiVar[1]][i],dxdyPhi(vect_[neiVar[1]], nbNeighbor[neiVar[1]], vect_[neiVar[0]], nbNeighbor[neiVar[0]], a[i], dxdyf_, g_, dxg_, type_)))),
	                              (multiple(grDxPhi1,vf_la_[neiVar[1]][i]))),
	                              (multiple(gra_vf_la_[neiVar[1]][i],dxPhi(vect_[neiVar[0]], nbNeighbor[neiVar[0]], vect_[neiVar[1]], nbNeighbor[neiVar[1]], a[i], dxf_, g_, dxg_, type_))));
	    vector<double>().swap(grPhi); // Free memory
	    vector<double>().swap(grDxPhi1); // Free memory
	    vector<double>().swap(grDxPhi2); // Free memory
	    vector<double>().swap(grDxDyPhi); // Free memory    
	  }
	  vector<int>().swap(neiVar); // Free memory
	}

  // Messages from variables to functions
    for (int unsigned i = 0; i < nVar; i++)
      if (!isLeaf_v[i]) {
	    vector <int> neiFunc; // Neighbor functions of variable
	    for (int unsigned j = 0; j < nFunc; j++) 
		  if (binMat_[i][j] == 1) 		   
		    neiFunc.push_back(j);
	    for (unsigned j = 0; j < neiFunc.size(); j++) {
	      double prod = 1;
	      double sum = 0;
	      vector <double> gra_sum(nPara), gra_sum2(nPara);
	      for (unsigned k = 0; k < neiFunc.size(); k++)
	        if ((k != j)&&(fv_mu_[i][neiFunc[k]]!=0)) {
	          prod *= fv_mu_[i][neiFunc[k]];
 		      sum += fv_la_[i][neiFunc[k]] / fv_mu_[i][neiFunc[k]];
 		      if (out_ > 0){
                gra_sum = add(multiple(gra_fv_mu_[i][neiFunc[k]], 1/fv_mu_[i][neiFunc[k]]), gra_sum);
                gra_sum2 = add(multiple(add(multiple(gra_fv_la_[i][neiFunc[k]], fv_mu_[i][neiFunc[k]]),
                               multiple(gra_fv_mu_[i][neiFunc[k]],- fv_la_[i][neiFunc[k]])), pow(fv_mu_[i][neiFunc[k]], -2)), gra_sum2);
              }
		    }
		  // Update lamda and mu functions from variable to function
		  vf_mu_[i][neiFunc[j]] = prod;
		  vf_la_[i][neiFunc[j]] = prod * sum;
		  if (out_ > 0){
		    gra_vf_mu_[i][neiFunc[j]] = multiple(gra_sum, prod);
		    gra_vf_la_[i][neiFunc[j]] = add(multiple(gra_vf_mu_[i][neiFunc[j]], sum),multiple(gra_sum2, prod));
		    vector<double>().swap(gra_sum); // Free memory
		    vector<double>().swap(gra_sum2); // free memory
	      }
		}
		vector<int>().swap(neiFunc); // Free memory
	  }
    nIteration_ --;
    if (nIteration_ < 1)
      isStop = true;
  }
  vector<double>().swap(theta_); // Free memory
  vector<double>().swap(vect_); // Free memory
  vector<int>().swap(nbNeighbor); // Free memory
  vector< vector<double> >().swap(a); // Free memory
  
  vector <int> neiRoot; // Neighbor functions of root
  for (unsigned i = 0; i < binMat_[root_].size(); i++)
    if (binMat_[root_][i] == 1)
      neiRoot.push_back(i);
  vector< vector<int> >().swap(binMat_); // Free memory
  
  // Initialisation values for results
  double u = 1;
  double z = 0;
  vector <double> gra_u(nPara), gra_z(nPara), gra_p(nPara);
    
  // Calculate P and gradient of P
  for (unsigned i = 0; i < neiRoot.size(); i++) {
    if (fv_mu_[root_][neiRoot[i]]!=0){
      u *= fv_mu_[root_][neiRoot[i]];
      z += fv_la_[root_][neiRoot[i]] / fv_mu_[root_][neiRoot[i]];
      if (out_ > 0){
        gra_u = add(multiple(gra_fv_mu_[root_][neiRoot[i]], 1/fv_mu_[root_][neiRoot[i]]), gra_u);
        gra_z = add(multiple(add(multiple(gra_fv_la_[root_][neiRoot[i]], fv_mu_[root_][neiRoot[i]]),
                    multiple(gra_fv_mu_[root_][neiRoot[i]], -fv_la_[root_][neiRoot[i]])), pow(fv_mu_[root_][neiRoot[i]] , -2)), gra_z);
      }
    }
  }
  vector< vector< vector<double> > >().swap(gra_fv_mu_); // Free memory
  vector< vector< vector<double> > >().swap(gra_fv_la_); // Free memory
  vector< vector< vector<double> > >().swap(gra_vf_mu_); // Free memory
  vector< vector< vector<double> > >().swap(gra_vf_la_); // free memory
  vector< vector<double> >().swap(fv_mu_); // Free memory
  vector< vector<double> >().swap(fv_la_); // Free memory
  vector< vector<double> >().swap(vf_mu_); // Free memory
  vector< vector<double> >().swap(vf_la_); // Free memory
  vector<int>().swap(neiRoot); // Free memory
  double p = u * z; // Compute density
  ret[0] = p;
  if (out_ > 0) {
    gra_u = multiple(gra_u, u);
    gra_p = add(multiple(gra_z, u),multiple(gra_u, z)); // Compute gradient
    for (int i = 1; i < ret.size(); i++)
      ret[i] = gra_p[i-1];
    vector<double>().swap(gra_u); // Free memory
    vector<double>().swap(gra_z); // Free memory
    vector<double>().swap(gra_p); // Free memory
    vector<bool>().swap(isLeaf_v); // Free memory
  }
  return ret;
}
