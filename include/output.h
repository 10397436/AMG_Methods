/**
* @file   output.h
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#ifndef OUTPUT_H_INCLUDED
#define OUTPUT_H_INCLUDED

#include "parameter_setup.h"
#include "parameter_cycle.h"
#include "parameter_method.h"

/** @class output
* @brief This class contains the printing and saving tools.
*
*
*/

class output
{
public:

/**
* @brief Constructor (defaulted)
*
*/

output()=default;

/**
* @brief Constructor
* @param[in] testname: name of output file
* @param[in] inputA: name of input file for matrix
* @param[in] inputf: name of input file for vector
* @param[in] fem: type of finite element discretization (CG conforming Galerkin, DG discontinuous Galerkin)
* @param[in] method: type of AMG method (AMG as stand-alone, PCG as preconditioner for conjugate gradient)
* @param[in] iter: number of iterations to achieve convergence
* @param[in] rho: convergence factor
* @param[in] flag: flag associated with convergence (0 convergence, 1 otherwise)
* @param[in] parameter_setup: parameters of setup
* @param[in] parameter_cycle: parameters of cycle
* @param[in] parameter_method: parameters of method
*
*/

output(const string& testname,const string& inputA, const string& inputf, const string& fem, const string& method, const int& iter,const Real& rho, const bool& flag, const parameter_setup& ps, const parameter_cycle& pc, const parameter_method& pm);

/**
* @brief Destructor (defaulted)
*
*/

~output(){};

/**
* @brief Print results on screen
*
*/

void print_on_screen() const;

/**
* @brief Print results on file
* @param[in] directory: file location path
*
*/

void print_on_file(const string& directory) const;

private:
string _testname; /**< @brief name of output file */ 
string _inputA; /**< @brief name of input file for matrix */
string _inputf; /**< @brief name of input file for vector */
string _fem; /**< @brief type of finite element discretization (CG conforming Galerkin, DG discontinuous Galerkin) */
string _method; /**< @brief type of AMG method (AMG as stand-alone, PCG as preconditioner for conjugate gradient) */
int _iter; /**< @brief number of iterations to achieve convergence */
bool _flag; /**< @brief flag associated with convergence (0 convergence, 1 otherwise) */
Real _rho; /**< @brief convergence factor */
parameter_setup _ps; /**< @brief parameters of setup */
parameter_cycle _pc; /**< @brief parameters of cycle */
parameter_method _pm; /**< @brief parameters of method */
};

#endif // OUTPUT_INCLUDED
