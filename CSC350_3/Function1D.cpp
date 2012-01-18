/*
 *  Function1D.cpp
 *  350Assignment1
 *
 *  Created by Los on 4/6/11.
 *  Copyright 2011 Universtity of Rhode Island. All rights reserved.
 *
 */
#include "Function1D.h"

namespace csc350Lib_calc_base {
	
	//Richardson Extrapolation 
	float Function1D::dfunc(float x)  
	{	
		float h = 1;
		float D[10][10];
		for (int i =0; i < 10; i++) {
			for (int j = 0; j<= i; j++) {
				if (j==0) {
					D[i][j] = (func(x+h)-func(x))/h;
				}
				else {
					D[i][j] = D[i][j-1] + (D[i][j-1] - D[i-1][j-1])/(pow(4, j) -1);
				}
				
			}
			h = h/2;
		}
		
		return D[9][9];	
		
	}//End of dfunc
	
}//end calc_base namespace