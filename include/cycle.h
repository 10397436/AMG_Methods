/**
* @file   cycle.h
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#ifndef CYCLE_H_INCLUDED
#define CYCLE_H_INCLUDED

#include "parameter_cycle.h"
#include "setup.h"

/** @class cycle
* @brief This class defines one iteration of mu-cycle.
*
*
*/

class cycle
{
public:

/**
* @brief Constructor (defaulted)
*
*/

cycle()=default;

/**
* @brief Constructor
* @param[in] S: setup containing coarser matrices and interpolation operators
* @param[in] f: right-hand side on finest level
* @param[in] p: parameters of cycle
*
*/

cycle(const setup& S, const Vec& f, const parameter_cycle& p); 

/**
* @brief Destructor (defaulted)
*
*/

~cycle(){};

/**
* @brief Reading setup S 
* @param[out] S: setup containing coarser matrices and interpolation operators
*
*/

inline const setup& get_S() const
{
return _S;
}

/**
* @brief Reading vector u 
* @param[in] n: current level of matrix/vector
* @param[out] u[n]: vector u at current level
*
*/

inline const Vec& get_u(const size_t& n)
{
if (n >= _pc.get_nlevel())  
{
	throw out_of_range("Index out of range.");
}
return _u[n];
}

/**
* @brief Writing vector u 
* @param[in] n: current level of matrix/vector
* @param[in] v: vector to be copied
* @param[out] u[n]: assigned vector u at current level
*
*/

inline void set_u(const size_t& n,const Vec& v)
{
if (n >= _pc.get_nlevel()) 
{
	throw out_of_range("Index out of range.");
}
_u[n]=v;
}

/**
* @brief Reading vector f
* @param[in] n: current level of matrix/vector
* @param[out] f[n]: vector f at current level
*
*/

inline const Vec& get_f(const size_t& n)
{
if (n >= _pc.get_nlevel()) 
{
throw out_of_range("Index out of range.");
}
return _f[n];
}

/**
* @brief Writing vector f
* @param[in] n: current level of matrix/vector
* @param[in] v: vector to be copied
* @param[out] f[n]: assigned vector f at current level
*
*/

inline void set_f(const size_t& n,const Vec& v)
{
if (n >= _pc.get_nlevel())  
{
throw out_of_range("Index out of range.");
}
_f[n]=v;
}

/**
* @brief Cycle iteration
* @param[in] lev: current level of coarser matrices, 0 is for finest level
*
*/

void Cycle(int lev);

private:

/**
* @brief Gauss-Seidel method
* @param[in] u: initial solution guess
* @param[in] j: current level of matrix/vector
* @param[in] f: right-hand side on j level
* @param[in] maxit: maximum number of smoothing iterations
*
*/

void GS(Vec& u,const Vec& f, const int& j, const int& maxit);

setup _S; /**< @brief setup containing coarser matrices and interpolation operators */
vector<Vec> _f; /**< @brief vector of right-hand side on all levels */
vector<Vec> _u; /**< @brief vector of solution on all levels */
parameter_cycle _pc; /**< @brief parameters of cycle */
};


#endif // CYCLE_H_INCLUDED
