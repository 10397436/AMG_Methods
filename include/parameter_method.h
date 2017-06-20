/**
* @file   parameter_method.h
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#ifndef PARAMETER_METHOD_H_INCLUDED
#define PARAMETER_METHOD_H_INCLUDED

#include "common.h"

/** @class parameter_method
* @brief This class contains method parameters.
*
*
*/

class parameter_method
{
public:

/**
* @brief Constructor (defaulted)
*
*/

parameter_method()=default;

/**
* @brief Constructor
* @param[in] tol: tolerance
* @param[in] maxiter: number of maximum iterations
*
*/

parameter_method(const Real& tol, const int& maxiter);

/**
* @brief Destructor (defaulted)
*
*/

~parameter_method(){};

/**
* @brief Reading parameter maxiter
* @param[out] maxiter: number of maximum iterations 
*
*/

inline const int& get_maxiter() const
{
return _maxiter;
}

/**
* @brief Reading parameter tol 
* @param[out] tol: tolerance
*
*/

inline const Real& get_tol() const
{
return _tol;
}

private:
Real _tol; /**< @brief tolerance */
int _maxiter; /**< @brief number of maximum iterations */
};

#endif // PARAMETER_METHOD_H_INCLUDED
