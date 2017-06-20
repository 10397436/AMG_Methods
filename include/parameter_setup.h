/**
* @file   parameter_setup.h
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#ifndef PARAMETER_SETUP_H_INCLUDED
#define PARAMETER_SETUP_H_INCLUDED

#include "common.h"

/** @class parameter_setup
* @brief This class contains setup parameters.
*
*
*/

class parameter_setup
{
public:

/**
* @brief Constructor (defaulted)
*
*/

parameter_setup()=default;

/**
* @brief Constructor
* @param[in] nmatrix: number of coarser matrices
* @param[in] theta: strong connection threshold
*
*/

parameter_setup(const int& nmatrix,const Real& theta);

/**
* @brief Destructor (defaulted)
*
*/

~parameter_setup(){};

/**
* @brief Reading parameter nmatrix
* @param[out] nmatrix: number of coarser matrices 
*
*/

inline const int& get_nmatrix() const
{
return _nmatrix;
}

/**
* @brief Reading parameter theta 
* @param[out] theta: strong connection threshold
*
*/

inline const Real& get_theta() const
{
return _theta;
}

private:
int _nmatrix; /**< @brief number of coarser matrices */
Real _theta; /**< @brief strong connection threshold */
};

#endif // PARAMETER_SETUP_H_INCLUDED
