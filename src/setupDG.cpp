/**
* @file   setupDG.cpp
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#include "setupDG.h"

setupDG::setupDG(const SpMat& A,const parameter_setup& p)
{
_ps=p;
_A.reserve(_ps.get_nmatrix());
_I.reserve(_ps.get_nmatrix()-1);
_A.push_back(A);
DG_setup();
}

void setupDG::DG_setup()
{
vector<sets> B;
aggregation_DG(B);
SpMat I(_A[0].rows(),B.size());
unsmoothed_interpolation(I,B);
GS_orth_interpolation(I);
smoothed_interpolation(I);
_I.push_back(I);
_A.push_back(I.transpose()*_A[0]*I);

for(int k=1;k<_ps.get_nmatrix();k++)
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

void setupDG::aggregation_DG(vector<sets>& B) 
{
sets currentB;
int R=_A[0].rows();
sets delset;
vector<Real> maxrow;
vector<int> pos;
maxrow_pos(_A[0],pos); //initialization
currentB.addElement(0);
currentB.addElement(pos[0]);
B.push_back(currentB);
currentB.clear_set();

for(int i=1;i<R;i++)
{
	int N=find_set(B,i);
	int M=find_set(B,pos[i]);
	if(_A[0].coeff(i,pos[i])>0)
	{
		if(N==-1) //new candidate singleton aggregate
		{
			if(delset.isEmpty())
			{
				currentB.addElement(i);
				B.push_back(currentB);
				currentB.clear_set();
			}
			else
			{
				delset.sort_set();
				B[delset[0]].addElement(i);
				delset.deleteElement(delset[0]);
			}
		}
	}
	else
	{
		if(N==-1 && M==-1)  //new candidate pair aggregate
		{
			if(delset.isEmpty())
			{
				currentB.addElement(i);
				currentB.addElement(pos[i]);
				B.push_back(currentB);
				currentB.clear_set();
			}
			else
			{
				delset.sort_set();
				B[delset[0]].addElement(i);
				B[delset[0]].addElement(pos[i]);
				delset.deleteElement(delset[0]);
			}
		}
		else //enlarging singleton or pair aggregates
		{
			if(N==-1 && M>=0) 
			{
				B[M].addElement(i);
			}
			else if(N>=0 && M==-1) 
			{
				B[N].addElement(pos[i]);
			}
			else if(N>=0 && M>=0 && M!=N) 
			{
				B[min(N,M)]=sets::union_set(B[N],B[M]);
				B[max(N,M)].clear_set();
				delset.addElement(max(N,M));
			}
		}
	}
}
}

void setupDG::unsmoothed_interpolation(SpMat& I, vector<sets>& B)
{
vector<Trip> elements;
for(size_t i=0;i<B.size();i++)
{
	for(int j=0;j<B[i].cardinality();j++)
		elements.push_back(Trip(B[i][j],i,1));
}

I.setFromTriplets(elements.begin(), elements.end()); //tentative interpolation operator
}

void setupDG::GS_orth_interpolation(SpMat& I)
{
//orthogonalization is not needed because columns of I are orhogonal by construction

//normalization
for(int i=0;i<I.cols();i++)
	I.col(i)=I.col(i)/(I.col(i)).norm();
}

void setupDG::smoothed_interpolation(SpMat& I)
{
Vec d=(_A[0].diagonal()).cwiseInverse();
SpMat D(_A[0].rows(),_A[0].cols());
D=d.asDiagonal();
SpMat Id(_A[0].rows(),_A[0].cols());
Id.setIdentity();
Real w=2./3;
I=(Id-w*D*_A[0])*I; //smoothing step
}

int setupDG::find_set(vector<sets>& B,const int& k)
{
int N=B.size();
for(int i=0;i<N;i++)
{
	if(B[i].isMember(k))
		return i;
}

return -1;
}

void setupDG::maxrow_pos(const SpMat& A, vector<int>& pos)
{
int dim=A.rows();
int ind;
Vec r;
Real val;
for(int i=0;i<dim;i++)
{
	r=(A.row(i)).cwiseAbs();
	r[i]=0;
	val=r.maxCoeff(&ind);
		if(val==0)
		{
			throw runtime_error("Possibly non DG matrix, found isolated point.");
		} 
	pos.push_back(ind);
}
}
