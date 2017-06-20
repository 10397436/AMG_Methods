/**
* @file   method.cpp
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#include  "method.h"

method::method(cycle& C,const parameter_method& p)
{
_C=C;
_pm=p;
_flag=0;
_iter=0;
_rho=0;
_solution=Vec::Zero((_C.get_f(0)).size());
}

void method::AMGCycle()
{
_iter=0;
_flag=0;
int init_lev(0);
Vec r=_C.get_f(init_lev)-_C.get_S().get_A(init_lev)*_C.get_u(init_lev);
Real r0=r.norm();
Real fnorm=(_C.get_f(init_lev)).norm();

if(fnorm==0)
	fnorm=1;

Real err=r0/fnorm;

while (err>_pm.get_tol() && _iter<_pm.get_maxiter()) //call to mu-cycle until convergence
{
	_iter++;
	_C.Cycle(init_lev);
	r=_C.get_f(init_lev)-_C.get_S().get_A(init_lev)*_C.get_u(init_lev);
	err=r.norm()/fnorm;
}

_solution=_C.get_u(init_lev);
Real rN=r.norm();
_rho=exp(log(rN/r0)/_iter);

if(err>_pm.get_tol()) //check if method is convergent
	_flag=1;
}

void method::PCGCycle()
{
_iter=0;
_flag=0;
_solution=Vec::Zero((_C.get_f(0)).size());
int init_lev(0);
Vec r=_C.get_f(init_lev)-_C.get_S().get_A(init_lev)*_solution;
Real r0=r.norm();
Real fnorm=(_C.get_f(init_lev)).norm();

if(fnorm==0)
	fnorm=1;

_C.set_u(init_lev,Vec::Zero((_C.get_f(init_lev)).size()));
_C.Cycle(init_lev); 
Real err=(_C.get_u(init_lev)).norm()/fnorm;
Real alpha=0,beta=0,csi=1e10,csiold=0;
Vec p,q;

//definition of preconditioned conjugate gradient method
while (err>_pm.get_tol() && _iter<_pm.get_maxiter()) //call to PCG until convergence
{
	_iter++;
	csi=r.transpose()*_C.get_u(init_lev);

	if(_iter>1)
	{
		beta=csi/csiold;
		p=_C.get_u(init_lev)+beta*p;
	}
	else
	{
		p=_C.get_u(init_lev);
	}

	q=_C.get_S().get_A(init_lev)*p;
	alpha=csi/(p.transpose()*q);
	_solution=_solution+alpha*p;
	r=r-alpha*q;
	_C.set_f(init_lev,r);
	_C.set_u(init_lev,Vec::Zero((_C.get_f(init_lev)).size()));
	_C.Cycle(init_lev); //AMG as preconditioner on residual equation
	err=(_C.get_u(init_lev)).norm()/fnorm;
	csiold=csi;
}

Real rN=r.norm();
_rho=exp(log(rN/r0)/_iter);

if(err>_pm.get_tol()) //check if method is convergent
	_flag=1;
}
