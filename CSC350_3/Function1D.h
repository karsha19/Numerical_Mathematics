/*
 *  Function1D.h
 *  350Assignment1
 *
 *  Created by Los on 4/6/11.
 *  Copyright 2011 Universtity of Rhode Island. All rights reserved.
 *
 */
#ifndef FUNCTION1D_H
#define FUNCTION1D_H
#include <math.h>

namespace csc350Lib_calc_base {
	
	/*This class is virtual. Never instanciate it directly.
	Instead, instanciate one of its hard-coded subclasses
	These subclasses will not belong to any package. They belong to the "application side". 
	Such classes would be implemented by a user of your library, who would want (in the future) 
	to use it to solve a nonlinear for his/her function. These functions will be hard-coded. 
	What the subclass does is simply to overload
	 always float func(float x)
	 possibly float dfunc(float x),[If you know a formula for the function's derivative]
	 always bool isExactDerivativeDefined(void)*/
	class Function1D {
	public:
		//Default constructor
		Function1D(void){};  
		
		//Always declare virtual the destructor of a class that contains a pure virtual function, 
		//even if you actually provide an implementation for the destructor. 
		//This takes care of potential problems related to polymorphism.
		virtual ~Function1D(void){}; 
		
		// This pure virtual function will be implemented by the subclasses.
		virtual float func(float x) = 0;    
		
		//This method should return the approximate value of your function's first derivative at x. 
		//This value will be computed using the Richardson extrapolation algorithm. 
		//Subclasses may overload this method, using a hard-coded expression for the exact derivative.
		virtual float dfunc(float x); 
		
		//This method, to be overloaded (and hard-coded) by subclasses, returns true if the subclass 
		//implements the exact derivatives, and false if the value returned by dfunc is computed 
		//using the Richardson extrapolation.
		virtual bool  isExactDerivativeDefined(void) = 0;
		
	};// end Function1D class
	
	
}//end calc_base namespace

#endif