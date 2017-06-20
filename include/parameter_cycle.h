/**
* @file   parameter_cycle.h
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#ifndef PARAMETER_CYCLE_H_INCLUDED
#define PARAMETER_CYCLE_H_INCLUDED

#include "common.h"

/** @class parameter_cycle
* @brief This class contains cycle parameters.
*
*
*/

class parameter_cycle
{
public:

/**
* @brief Constructor (defaulted)
*
*/

parameter_cycle()=default;

/**
* @brief Constructor
* @param[in] nlevel: number of coarser levels
* @param[in] nu1: number of pre-smoothing iterations
* @param[in] nu2: number of post-smoothing iterations
* @param[in] mu: flag to decide type of cycle: mu=1 V-cycle, mu=2 W-cycle
*
*/

parameter_cycle(const int& nlevel,const int& nu1,const int& nu2,const int& gamma);

/**
* @brief Destructor (defaulted)
*
*/

~parameter_cycle(){};

/**
* @brief Reading parameter nlevel 
* @param[out] nlevel: number of coarser levels
*
*/

inline const int& get_nlevel() const
{
return _nlevel;
}

/**
* @brief Reading parameter nu1
* @param[out] nu1: number of pre-smoothing iterations
*
*/

inline const int& get_nu1() const
{
return _nu1;
}

/**
* @brief Reading parameter nu2
* @param[out] nu2: number of post-smoothing iterations
*
*/

inline const int& get_nu2() const
{
return _nu2;
}

/**
* @brief Reading parameter mu 
* @param[out] mu: flag to decide type of cycle: mu=1 V-cycle, mu=2 W-cycle
*
*/

inline const int& get_mu() const
{
return _mu;
}

private:
int _nlevel; /**< @brief number of coarser levels */
int _nu1; /**< @brief number of pre-smoothing iterations */
int _nu2; /**< @brief number of post-smoothing iterations */
int _mu; /**< @brief flag to decide type of cycle: mu=1 V-cycle, mu=2 W-cycle */
};

#endif // PARAMETER_H_INCLUDED
