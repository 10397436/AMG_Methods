/**
* @file   setupDG.h
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#ifndef SETUPDG_H_INCLUDED
#define SETUPDG_H_INCLUDED

#include "sets.h"
#include "setup.h"
#include "parameter_setup.h"

/** @class setupDG
* @brief Class inherited from setup. This class defines the construction of coarser matrices and interpolation operators for matrices stemming from discontinous Galerkin discretization extending the algorithm for conforming Galerkin matrices.
*
*
*/

class setupDG: public setup
{
public:

/**
* @brief Constructor (defaulted)
*
*/

setupDG()=default;

/**
* @brief Constructor
* @param[in] A: input matrix defined on finest level
* @param[in] p: parameters of setup
*
*/

setupDG(const SpMat& A,const parameter_setup& p);

/**
* @brief Destructor (defaulted)
*
*/

~setupDG(){};

private:

/**
* @brief Aggregation
* @param[in] B: initialization of vector of aggregate sets (it will be built in the method)
*
*/

void aggregation_DG(vector<sets>& B);
/**
* @brief Unsmoothed interpolation formula
* @param[in] I: initialization of interpolation operator (it will be built in the method)
* @param[in] B: vector of aggregate sets
*
*/
void unsmoothed_interpolation(SpMat& I, vector<sets>& B);

/**
* @brief Gram-Schmidt orthonormalization applied to the interpolation formula
* @param[in] I: interpolation operator
*
*/

void GS_orth_interpolation(SpMat& I);

/**
* @brief Smoothing step applied to the interpolation formula
* @param[in] I: interpolation operator
*
*/

void smoothed_interpolation(SpMat& I);

/**
* @brief Construnction of coarser matrices and interpolation operators for matrix stemming from discontinuous Galerkin discretization 
*
*/

void DG_setup();

/**
* @brief Find the aggregate set containing a given value 
* @param[in] B: vector of aggregate sets
* @param[in] k: value to be found
* @param[out] i: index of set containing k, if it is not found then i=-1
*
*/

int find_set(vector<sets>& B,const int& k);

/**
* @brief Utility: 
* @param[in] A: input matrix defined on finest level
* @param[in] pos: initialization of vector containing position of all maximum values of all matrix rows except for diagonal values (it will be built in the method)
*
*/

void maxrow_pos(const SpMat& A, vector<int>& pos);
};

#endif // SETUPDG_H_INCLUDED
