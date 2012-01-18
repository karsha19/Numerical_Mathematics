/*
 *  PolyFunction1D.h
 *  350Ass1
 *
 *  Created by Los on 4/12/11.
 *  Copyright 2011 Universtity of Rhode Island. All rights reserved.
 *
 */
#ifndef POLYFUNCTION1D_H
#define POLYFUNCTION1D_H
#include "Function1D.h"

//This class simply implements polynomial functions. 
//A polynomial function is entirely defined by its list of coefficients, 
//which will be passed to the class constructor
class PolyFunction1D : public csc350Lib_calc_base::Function1D{
public: 
	float func(float);
	float trigonometric_func(float, float ,float);
	float dfunc(float);
	bool isExactDerivativeDefined(void);
	PolyFunction1D(int nbCoeffs, const float *coefficients);
	PolyFunction1D(){};
	~PolyFunction1D(){};
private:
	int numberOfCoeffs_;
	float coeff_[5];
};

#endif