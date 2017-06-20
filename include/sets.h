/**
* @file   sets.h
* @author Laura Melas <laura.melas@mail.polimi.it>
* @date   2017
*
* This file is part of project "AMG Methods".
*
* @brief AMG methods for conforming and discontinuous Galerkin finite element discretizations of the Poisson problem.
*
*/

#ifndef SETS_INCLUDED
#define SETS_INCLUDED

#include "common.h"

/** @class sets
* @brief This class performs some properties and utilities of mathematical sets.
*
*
*/

class sets
{
public:

/**
* @brief Constructor (defaulted)
*
*/

sets()=default;

/**
* @brief Constructor
* @param[in] dim: cardinality of the set
*
*/

sets(const size_t& dim);

/**
* @brief Copy constructor
* @param[in] A: set to be copied 
*
*/

sets(const sets& A);

/**
* @brief Destructor (defaulted)
*
*/

~sets(){}

/**
* @brief Definition of operator [], writing version
* @param[in] n: access position to an element of the set 
* @param[out] set[n]: write element in the choosen position of the set 
*
*/

inline int& operator[](const size_t& n){
if (n >= _set.size()) 
{
	throw out_of_range("Index out of range.");
}
return _set[n];
}

/**
* @brief Definition of operator [], reading version
* @param[in] n: access position to an element of the set 
* @param[out] set[n]: read element in the choosen position of the set 
*
*/

inline const int& operator[](const size_t& n) const{
if (n >= _set.size()) 
{
	throw out_of_range("Index out of range.");
}
return _set[n];
}

/**
* @brief Add an element in the set
* @param[in] s: element to be added
*
*/

void addElement(const int& s);

/**
* @brief Delete an element in the set
* @param[in] s: element to be deleted
*
*/

void deleteElement(const int& s);

/**
* @brief Check if an element is in the set
* @param[in] s: element to be found
* @param[out] 0,1    : 1 if s is in the set, 0 otherwise
*
*/

bool isMember(const int& s);

/**
* @brief Check if a set is empty
* @param[out] 0,1    : 1 if the set is empty, 0 otherwise
*
*/

bool isEmpty();

/**
* @brief Find the position of an element in the set
* @param[in] s: element to be found
* @param[out] d : position of the element
*
*/

int find_pos_set(const int& s);

/**
* @brief Cardinality of the set
*
*/

int cardinality();

/**
* @brief Union between two sets
* @param[out] U    : union set between A and B
* @param[in]  A,B        : two sets
*
*/

static sets union_set(sets& A,sets& B);

/**
* @brief Difference between two sets 
* @param[out] D    : difference set between A and B (D=A-B)
* @param[in]  A,B        : two sets
*
*/

static sets diff_set(sets& A,sets& B);

/**
* @brief Intersection between two sets
* @param[out] I    : intersection set between A and B
* @param[in]  A,B        : two sets
*
*/

static sets inter_set(sets& A,sets& B);

/**
* @brief Reorder the set
*
*/

void sort_set();

/**
* @brief Delete all element of the set
*
*/

void clear_set();

private:
vector<int> _set; /**< @brief definition of set */
};

#endif // SETS_INCLUDED
