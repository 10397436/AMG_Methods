/**
* @file   sets.cpp
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#include "sets.h"

sets::sets(const size_t& dim):_set(dim) {}

sets::sets(const sets& A) : _set(A._set) {}

bool sets::isMember(const int& s)
{
return binary_search(_set.begin(),_set.end(),s);
}

int sets::cardinality()
{
return _set.size();
}

void sets::addElement(const int& s)
{
_set.push_back(s);
}

int sets::find_pos_set(const int& s)
{
auto it=find(_set.begin(),_set.end(),s);
int d=distance(_set.begin(),it);
return d;
}

void sets::deleteElement(const int& s)
{
_set.erase(find(_set.begin(),_set.end(),s));
}

void sets::sort_set()
{
sort(_set.begin(),_set.end());
}

sets sets::union_set(sets& A,sets& B)
{
sets U(A.cardinality()+B.cardinality());
A.sort_set();
B.sort_set();
auto it=set_union((A._set).begin(),(A._set).end(),(B._set).begin(),(B._set).end(),(U._set).begin());
(U._set).resize(it-(U._set).begin());;
return U;
}

sets sets::diff_set(sets& A,sets& B)
{
sets D(A.cardinality()+B.cardinality());
A.sort_set();
B.sort_set();
auto it=set_difference((A._set).begin(),(A._set).end(),(B._set).begin(),(B._set).end(),(D._set).begin());
(D._set).resize(it-(D._set).begin());
return D;
}

sets sets::inter_set(sets& A,sets& B)
{
sets I(A.cardinality()+B.cardinality());
A.sort_set();
B.sort_set();
auto it=set_intersection((A._set).begin(),(A._set).end(),(B._set).begin(),(B._set).end(),(I._set).begin());
(I._set).resize(it-(I._set).begin());
return I;
}

bool sets::isEmpty()
{
return _set.empty();
}

void sets::clear_set()
{
_set.clear();
}

