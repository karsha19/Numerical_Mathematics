/*
 *  PolyFunction1D.cpp
 *  350Ass1
 *
 *  Created by Los on 4/12/11.
 *  Copyright 2011 Universtity of Rhode Island. All rights reserved.
 *
 */
#define PI 3.14159
#include <iostream>
#include "PolyFunction1D.h"
using namespace std;

//Given the number of coefficients and the coefficients model a polynomial
PolyFunction1D::PolyFunction1D(int nbCoeffs, const float * coefficients)
{ 
	numberOfCoeffs_ = nbCoeffs;
	for (int i =0; i < nbCoeffs; i++) {
		coeff_[i] = *(coefficients + i);
	}
}
//Return true because the exact derivative is defined
bool PolyFunction1D::isExactDerivativeDefined(void)
{
	return true;
}



float PolyFunction1D::trigonometric_func(float x, float L, float n)
{
	
	float temp = 0;
	
	for (int i = 0; i< numberOfCoeffs_; i++) {
		if (i == 0) {
			temp = coeff_[i];
		}
		else {
			for (int columns = 1; columns <= n; columns++) {
				temp += coeff_[i] + cos((2 * columns * PI * x)/ L );
				
			}
			for (int columns = 1; columns <= n; columns++) {
				temp += coeff_[i] + sin((2 * columns * PI * x)/ L );
			}
			
		}

	}
	return temp;

}

//Find the value of the polynomial at a given x
float PolyFunction1D::func(float x)
{
	float temp=0;
	for (int i=0; i<numberOfCoeffs_; i++) {
		temp = temp + coeff_[i] * pow(x, i);
	}
	return temp;
}


//Find the derivative of polynomial at a given point (x) using power method
float PolyFunction1D::dfunc(float x)
{
	float temp=0;
	for (int i =0; i<numberOfCoeffs_; i++) {
		temp = temp + i * coeff_[i] * pow(x, i-1);
	}
	return temp;
}
