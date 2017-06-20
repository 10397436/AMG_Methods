/**
* @file   parameter_setup.cpp
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#include "parameter_setup.h"

parameter_setup::parameter_setup(const int& nmatrix,const Real& theta)
{
_nmatrix=nmatrix;
_theta=theta;
}


