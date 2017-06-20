/**
* @file   parameter_cycle.cpp
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#include "parameter_cycle.h"

parameter_cycle::parameter_cycle(const int& nlevel,const int& nu1,const int& nu2,const int& mu)
{
_nlevel=nlevel;
_nu1=nu1;
_nu2=nu2;
_mu=mu;
}

