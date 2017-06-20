/**
* @file   parameter_method.cpp
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#include "parameter_method.h"

parameter_method::parameter_method(const Real& tol, const int& maxiter)
{
_tol=tol;
_maxiter=maxiter;
}

