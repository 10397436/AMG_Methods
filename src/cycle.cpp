/**
* @file   cycle.cpp
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#include "cycle.h"

cycle::cycle(const setup& S, const Vec& f, const parameter_cycle& p)
{
_S=S;
_pc=p;
_u.resize(_pc.get_nlevel()+1);
_f.resize(_pc.get_nlevel()+1);
_u[0]=Vec::Zero(f.size());
_f[0]=f;
}

void cycle::GS(Vec& u,const Vec& f, const int& j, const int& maxit)
{
auto L=_S.get_A(j).triangularView<Lower>();
Vec r=f-_S.get_A(j)*u;
Vec z;
int iter(0);

while (iter<maxit)
{
	iter++;
	z=L.solve(r);
	u=u+z;
	r=r-_S.get_A(j)*z;
}
}

void cycle::Cycle(int lev)
{
GS(_u[lev],_f[lev],lev,_pc.get_nu1()); //pre-smoothing
if(lev==_pc.get_nlevel())
{
	SimplicialLLT<SpMat> solver;
	solver.analyzePattern(_S.get_A(lev));
	solver.factorize(_S.get_A(lev));
	_u[lev]=solver.solve(_f[lev]);  //direct solver
}
else
{
	_f[lev+1]=(_S.get_I(lev)).transpose()*(_f[lev]-_S.get_A(lev)*_u[lev]);
	_u[lev+1]=Vec::Zero(_f[lev+1].size());
	for(int c=0;c<_pc.get_mu();c++)  //recursive call of mu-cycle
	{
		Cycle(lev+1);
		if(lev+1==_pc.get_nlevel()) //this break is to avoid to solve twice the same linear system when mu=2
		break;
	}
	_u[lev]=_u[lev]+_S.get_I(lev)*_u[lev+1];
	GS(_u[lev],_f[lev],lev,_pc.get_nu2()); //post-smoothing
}
}

