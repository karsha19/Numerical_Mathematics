#ifndef	INTERPOLATION_H
#define INTERPOLATION_H
#include "PolyFunction1D.h"
#include "Matrix.h"
#include "ColVec.h"
#include <cmath>

class Interpolation{
public:
	//Default constructor
	Interpolation(void);
	
	//Destructor
	~Interpolation(void){};
	
	PolyFunction1D polynomial_case(int numberOfPoints, float pointsArray[][2] );
	PolyFunction1D trigonometric_case( int numberOfPoints, float pointsArray[][2], float, float );
	float lagrange_case(int numberOfPoints, float pointsArray[][2], float x );
private:
	

};

#endif