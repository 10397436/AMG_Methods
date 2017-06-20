/**
* @file   setup.cpp
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#include "setup.h"
#include "sets.h"

setup::setup(const SpMat& A,const parameter_setup& p)
{
_ps=p;
_A.reserve(_ps.get_nmatrix());
_I.reserve(_ps.get_nmatrix()-1);
_A.push_back(A);
CG_setup();
}

void setup::CG_setup()
{
for(int k=0;k<_ps.get_nmatrix();k++)
{
	size_t n=_A[k].rows();
	sets C,F;
	vector<sets> S(n),St(n),Ci(n),Ds(n),Dw(n);
	strong_influence_dependence(_A[k],S,St,Dw);
	colouring_scheme(S,St,C,F);
	coarse_strong_dependence(S,Ci,Ds,C);
	check_modify(C,F,Ci,Ds);
	SpMat I(n,C.cardinality());
	interpolation(_A[k],I,C,Ci,Ds,Dw);
	_I.push_back(I);
	_A.push_back(I.transpose()*_A[k]*I);
}
}

void setup::strong_influence_dependence(const SpMat& A, vector<sets>& S, vector<sets>& St, vector<sets>& Dw)
{
vector<Real> maxrow,maxcol;
minus_maxrow_maxcol(A,maxrow,maxcol);
for (int k=0; k < A.outerSize(); ++k)
{
	for (SpMat::InnerIterator it(A,k); it; ++it)
	{
		if(it.row()!=it.col())
		{
			if(abs(it.value())>=_ps.get_theta()*maxrow[it.row()]) //definition of strong dependence set
			{
				S[it.row()].addElement(it.col());
			}
			else
			{
				Dw[it.row()].addElement(it.col());
			}
			if(abs(it.value())>=_ps.get_theta()*maxcol[it.col()]) //definition of strong influence set
			{
				St[it.row()].addElement(it.col());
			}
		}
	}
}
}

void setup::minus_maxrow_maxcol(const SpMat& A,vector<Real>& maxrow, vector<Real>& maxcol)
{
int dim=A.rows();
Vec r;
for(int i=0;i<dim;i++)
{
	r=(A.row(i)).cwiseAbs();
	r[i]=0;
	maxrow.push_back(r.maxCoeff());
	r=(A.col(i)).cwiseAbs();
	r[i]=0;
	maxcol.push_back(r.maxCoeff());
}
}

void setup::colouring_scheme(vector<sets>& S, vector<sets>& St, sets& C, sets& F)
{
int n,m,p;
size_t I;
vector<int> lambda;
for(size_t i=0;i<St.size();i++)
{
	lambda.push_back(St[i].cardinality()); //measure lambda
}
sets newF,U,J;
while(*max_element(lambda.begin(),lambda.end())!=-1)
{
	I=distance(lambda.begin(),max_element(lambda.begin(),lambda.end())); //new C point
	C.addElement(I);
	newF=St[I]; 
	U=sets::union_set(C,F);
	newF=sets::diff_set(newF,U); //new F points
	n=newF.cardinality();
	for(int i=0;i<n;i++)
	{
		F.addElement(newF[i]);
	}
	U=sets::union_set(C,F);
	m=U.cardinality();
	for(int i=0;i<m;i++) //update lambda
	{
		lambda[U[i]]=-1;
	}
	for(int i=0;i<n;i++)
	{
		J=S[newF[i]];
		J=sets::diff_set(J,U);
		p=J.cardinality();
		for(int j=0;j<p;j++)
		{
			++lambda[J[j]];
		}
	}
}

C.sort_set();
F.sort_set();
}

void setup::coarse_strong_dependence(vector<sets>& S, vector<sets>& Ci, vector<sets>& Ds, sets C)
{
int n=S.size();
int m;
for(int i=0;i<n;i++)
{
	m=S[i].cardinality();
	for(int j=0;j<m;j++) //definition of coarse interpolatory set and strong non-interpolatory set
	{
		if(C.isMember(S[i][j]))
			Ci[i].addElement(S[i][j]);
		else
			Ds[i].addElement(S[i][j]);
	}
}
}

void setup::check_modify(sets& C, sets& F, vector<sets>& Ci, vector<sets>& Ds)
{
size_t n=F.cardinality();
size_t all=Ci.size();
size_t m;
bool b;

for(size_t i=0;i<n;i++)
{
	m=Ds[F[i]].cardinality();
	b=0;
	for(size_t j=0;j<m;j++)
	{
		if((sets::inter_set(Ci[F[i]],Ci[Ds[F[i]][j]])).isEmpty() && b==0) //check if heuristics are satisfied
		{
			b=1;
			C.addElement(F[i]);
			for(size_t k=0;k<all;k++)
			{
				if(Ds[k].isMember(F[i]))
				{
					Ci[k].addElement(F[i]);
					Ds[k].deleteElement(F[i]);
				}
			}
		}
	}
}

F=sets::diff_set(F,C);
C.sort_set();
F.sort_set();
}

void setup::interpolation(const SpMat& A, SpMat& I, sets& C, const vector<sets>& Ci, const vector<sets>& Ds, const vector<sets>& Dw)
{
Real g;
size_t N=A.rows();
vector<Trip> elements;
for(size_t i=0;i<N;i++)
{
	if(C.isMember(i))
	{
		size_t pos=C.find_pos_set(i);
		elements.push_back(Trip(i,pos,1));
	}
	else //computation of interpolation weigths
	{
		sets Cii=Ci[i],Dis=Ds[i],Diw=Dw[i];
		size_t c=Cii.cardinality(),s=Dis.cardinality(),w=Diw.cardinality();
		Real den=A.coeff(i,i);
		SpVec X(N),E(N),S(N);
		SpCount L(N);

		//compute denominator
		for(size_t n=0;n<w;n++)
		{
			Real sum(0),sumabs(0);
			vector<Real> Aeval=element_set(A,Cii,Diw[n]);
			for(auto it=Aeval.begin();it != Aeval.end(); ++it)
			{
				sum+=*it;
				sumabs+=abs(*it);
			}
			L.insert(Diw[n])=Aeval.size();
			X.insert(Diw[n])=-sum/sumabs;
			S.insert(Diw[n])=sumabs;
			if(L.coeff(Diw[n])==0)
			{
				den-=abs(A.coeff(i,Diw[n]));
			}
			else if(X.coeff(Diw[n])>=0.5 && A.coeff(i,Diw[n])<0)
			{
				den-=A.coeff(i,Diw[n]);
			}
		}

		for(size_t m=0;m<s;m++)
		{
			Real sum(0),sumabs(0);
			vector<Real> Aeval=element_set(A,Cii,Dis[m]);
			for(auto it=Aeval.begin();it != Aeval.end(); ++it)
			{
				sum+=*it;
				sumabs+=abs(*it);
			}
			L.insert(Dis[m])=Aeval.size();
			X.insert(Dis[m])=-sum/sumabs;
			S.insert(Dis[m])=sumabs;
			E.insert(Dis[m])=abs(A.coeff(Dis[m],i))*L.coeff(Dis[m])/sumabs;
			if(E.coeff(Dis[m])<0.75 && X.coeff(Dis[m])>=0.5 && A.coeff(i,Dis[m])<0)
			{
				den-=A.coeff(i,Dis[m]);
			}
			else if(E.coeff(Dis[m])>2 && X.coeff(Dis[m])>=0.5 && A.coeff(i,Dis[m])<0)
			{
				den+=0.5*A.coeff(i,Dis[m]);
			}
		}

		//compute numerator
		for(size_t j=0;j<c;j++)
		{
			Real num=A.coeff(i,Cii[j]);

			for(size_t n=0;n<w;n++)
			{
				if(L.coeff(Diw[n])>0 && X.coeff(Diw[n])>=0.5 && A.coeff(i,Diw[n])<0)
				{
					g=abs(A.coeff(Diw[n],Cii[j]))/S.coeff(Diw[n]);
					num+=2*g*A.coeff(i,Diw[n]);
				}
				else if(L.coeff(Diw[n])>0)
				{
					g=abs(A.coeff(Diw[n],Cii[j]))/S.coeff(Diw[n]);
					num+=g*A.coeff(i,Diw[n]);
				}
			}

			for(size_t m=0;m<s;m++)
			{
				if(E.coeff(Dis[m])<0.75 && X.coeff(Dis[m])>=0.5 && A.coeff(i,Dis[m])<0)
				{
					g=abs(A.coeff(Dis[m],Cii[j]))/S.coeff(Dis[m]);
					num+=2*g*A.coeff(i,Dis[m]);
				}
				else if(E.coeff(Dis[m])>2 && X.coeff(Dis[m])>=0.5 && A.coeff(i,Dis[m])<0)
				{
					g=abs(A.coeff(Dis[m],Cii[j]))/S.coeff(Dis[m]);
					num+=0.5*g*A.coeff(i,Dis[m]);
				}
				else
				{
					g=abs(A.coeff(Dis[m],Cii[j]))/S.coeff(Dis[m]);
					num+=g*A.coeff(i,Dis[m]);
				}
			}
			Real wij=-num/den;
			size_t pos=C.find_pos_set(Cii[j]);
			elements.push_back(Trip(i,pos,wij));	
		}
	}
}

I.setFromTriplets(elements.begin(), elements.end());
}

vector<Real> setup::element_set(const SpMat& A, sets& B, const int& c)
{
vector<Real> eval;
size_t N=B.cardinality();
	for(size_t i=0;i<N;i++)
	{
		if(A.coeff(c,B[i])!=0)
		{
			eval.push_back(A.coeff(c,B[i]));
		}
	}

return eval;
}

