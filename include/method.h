/**
* @file   method.h
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#ifndef METHOD_INCLUDED
#define METHOD_INCLUDED

#include  "cycle.h"
#include  "parameter_method.h"

/** @class method
* @brief This class defines AMG methods: AMG stand-alone or PCG preconditioned conjugate gradient.
*
*
*/

class method
{
public:

/**
* @brief Constructor (defaulted)
*
*/

method()=default;

/**
* @brief Constructor
* @param[in] C: definition of one iteration of mu-cycle
* @param[in] p: parameters of method
*
*/

method(cycle& C,const parameter_method& p);

/**
* @brief Destructor (defaulted)
*
*/

~method(){};

/**
* @brief AMG stand-alone method
*
*/

void AMGCycle();

/**
* @brief PCG method (AMG as preconditioner)
*
*/

void PCGCycle();

/**
* @brief Reading number of iterations to achieve convergence
* @param[out] iter: number of iterations to achieve convergence
*
*/

inline const int& get_iter() const
{
return _iter;
}

/**
* @brief Reading flag
* @param[out] flag: flag associated with convergence of method (0 convergence, 1 otherwise)
*
*/

inline const bool& get_flag() const
{
return _flag;
}

/**
* @brief Reading convergence factor
* @param[out] rho: convergence factor of method 
*
*/

inline const Real& get_rho() const
{
return _rho;
}

/**
* @brief Reading solution vector
* @param[out] solution: vector containing solution 
*
*/

inline const Vec& get_solution() const
{
return _solution;
}

private:
cycle _C; /**< @brief definition of one iteration of mu-cycle */
Vec _solution; /**< @brief vector containing solution */
int _iter; /**< @brief number of iterations to achieve convergence */
Real _rho; /**< @brief convergence factor of method */ 
bool _flag; /**< @brief flag associated with convergence of method (0 convergence, 1 otherwise) */
parameter_method _pm; /**< @brief parameters of method */
};


#endif // METHOD_INCLUDED
