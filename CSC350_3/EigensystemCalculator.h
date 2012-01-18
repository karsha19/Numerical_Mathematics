#ifndef EIGENSYSTEMCALCULATOR_H
#define EIGENSYSTEMCALCULATOR_H
#include "RealEigenpair.h"

class EigensystemCalculator{
public:
	//Find largest eigenpair of given Matrix
	const RealEigenpair *get_largest_eigenpair(const Matrix<float> *);
	
	//Find # of eigenpairs. Parameters: Matrix, # of eigenpairs to return
	const RealEigenpair* get_eigenpairs(const Matrix<float> *, int);
};


#endif EIGENSYSTEMCALCULATOR_H
