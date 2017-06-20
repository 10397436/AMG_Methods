/**
* @file   common.h
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/SparseExtra>
#include <iostream>
#include <fstream>
#include "stdlib.h"
#include <vector>
#include <algorithm> 
#include <iterator>
#include <cmath> 

using namespace Eigen;
using namespace std;
using Real=double; /**< @brief Typedef for real numbers. */
typedef SparseMatrix<Real> SpMat; /**< @brief Typedef for sparse real-valued matrices. */
typedef SparseVector<Real> SpVec; /**< @brief Typedef for sparse real-valued vectors. */
typedef SparseVector<int> SpCount; /**< @brief Typedef for sparse int-valued vectors. */
using Vec=Matrix<Real,Dynamic,1>;/**< @brief Typedef for real-valued vectors. */ 
typedef Triplet<Real> Trip; /**< @brief Typedef for triplet, it is used to build sparse real-valued matrices. */

#endif // COMMON_H_INCLUDED
