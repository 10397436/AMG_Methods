/**
* @file   setup.h
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#ifndef SETUP_H_INCLUDED
#define SETUP_H_INCLUDED

#include "sets.h"
#include "parameter_setup.h"

/** @class setup
* @brief This class defines the construction of coarser matrices and interpolation operators, in particular for matrices stemming from conforming Galerkin discretization.
*
*
*/

class setup
{
public:

/**
* @brief Constructor (defaulted)
*
*/

setup()=default;

/**
* @brief Constructor
* @param[in] A: input matrix defined on finest level
* @param[in] p: parameters of setup
*
*/

setup(const SpMat& A,const parameter_setup& p);

/**
* @brief Destructor (defaulted)
*
*/

~setup(){};

/**
* @brief Reading matrix A
* @param[in] n: current position of coarser matrix
* @param[out] A[n]: matrix A at current position
*
*/

inline const SpMat& get_A(const size_t& n) const
{
if (n >= _A.size()) {
	throw out_of_range("Index out of range.");
}
return _A[n];
}

/**
* @brief Reading interpolation operator I
* @param[in] n: current position of interpolation operator
* @param[out] I[n]: operator I at current position
*
*/

inline const SpMat& get_I(const size_t& n) const
{
if (n >= _I.size()) {
	throw out_of_range("Index out of range.");
}
return _I[n];
}

protected:

/**
* @brief Definition of strong connections
* @param[in] A: input matrix defined on finest level
* @param[in] S: initialization of vector of sets containing all strong dependence connections (it will be built in the method)
* @param[in] St: initialization of vector of sets containing all strong influence conncetions (it will be built in the method)
* @param[in] Dw: initialization of vector of sets containing all weak connections (it will be built in the method)
*
*/

void strong_influence_dependence(const SpMat& A, vector<sets>& S, vector<sets>& St, vector<sets>& Dw);

/**
* @brief First step of coarsening strategy: C/F splitting
* @param[in] S: vector of sets containing all strong dependence connections
* @param[in] St: vector of sets containing all strong influence conncetions
* @param[in] C: initialization of C-points (it will be built in the method)
* @param[in] F: initialization of F-points (it will be built in the method)
*
*/

void colouring_scheme(vector<sets>& S, vector<sets>& St, sets& C, sets& F);

/**
* @brief Definition of vectors of coarse-interpolatory sets and of strong non-interpolatory sets
* @param[in] S: vector of sets containing all strong dependence connections
* @param[in] Ci: initialization of vector of coarse interpolatory sets (it will be built in the method)
* @param[in] Ds: initialization of vector of strong non-interpolatory sets (it will be built in the method)
* @param[in] C: C-points of C/F-splitting
*
*/

void coarse_strong_dependence(vector<sets>& S, vector<sets>& Ci, vector<sets>& Ds, sets C);

/**
* @brief Second step of coarsening strategy: C/F splitting
* @param[in] C: C-points of C/F-splitting
* @param[in] F: F-points of C/F-splitting
* @param[in] Ci: vector of coarse interpolatory sets
* @param[in] Ds: vector of strong non-interpolatory sets
*
*/

void check_modify(sets& C, sets& F, vector<sets>& Ci, vector<sets>& Ds);

/**
* @brief Interpolation formula
* @param[in] A: input matrix defined on finest level
* @param[in] I: initialization of interpolation operator (it will be built in the method)
* @param[in] C: C-points of C/F-splitting
* @param[in] Ci: vector of coarse interpolatory sets
* @param[in] Ds: vector of strong non-interpolatory sets
* @param[in] Dw: vector of weak non-interpolatory sets
*
*/


void interpolation(const SpMat& A, SpMat& I, sets& C, const vector<sets>& Ci, const vector<sets>& Ds, const vector<sets>& Dw);

/**
* @brief Construnction of coarser matrices and interpolation operators for matrix stemming from conforming Galerkin discretization 
*
*/

void CG_setup();

/**
* @brief Utility: 
* @param[in] A: input matrix defined on finest level
* @param[in] B: indices set
* @param[in] c: index
* @param[out] Aeval: vector containing all values of A(B,c)
*
*/

vector<Real> element_set(const SpMat& A, sets& B, const int& c);

/**
* @brief Utility: 
* @param[in] A: input matrix defined on finest level
* @param[in] maxrow: initialization of vector containing maximum values of all matrix rows (it will be built in the method)
* @param[in] maxcol: initialization of vector containing maximum values of all matrix columns (it will be built in the method)
*
*/

void minus_maxrow_maxcol(const SpMat& A,vector<Real>& maxrow, vector<Real>& maxcol);

vector<SpMat> _A; /**< @brief vector containing coarser matrices */
vector<SpMat> _I; /**< @brief vector containing interpolation operators */
parameter_setup _ps; /**< @brief parameters of setup */
};

#endif // SETUP_H_INCLUDED
