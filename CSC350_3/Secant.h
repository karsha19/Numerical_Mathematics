/*
 *  Secant.h
 *  350Ass1
 *
 *  Created by Los on 4/13/11.
 *  Copyright 2011 Universtity of Rhode Island. All rights reserved.
 *
 */
#ifndef SECANT_H
#define SECANT_H
#include "SolutionNLE.h"
#include "NonlinearSolver.h"

class Secant : public csc350Lib_calculus_snle::NonlinearSolver {
	
public:
	Secant() : NonlinearSolver(){};
	//solve for one function
	csc350Lib_calculus_snle::SolutionNLE* solve(csc350Lib_calc_base::Function1D *f, float, float, float tol);
	
	//solve for two functions
	csc350Lib_calculus_snle::SolutionNLE* solve_for_2(csc350Lib_calc_base::Function1D *f,csc350Lib_calc_base::Function1D *f2, float x, float a, float tol);
	
	
	
};

#endif
