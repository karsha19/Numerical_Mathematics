/*
 *  Secant.cpp
 *  350Ass1
 *
 *  Created by Los on 4/13/11.
 *  Copyright 2011 Universtity of Rhode Island. All rights reserved.
 *
 */
#include "Secant.h"
#include <iostream>
using namespace std;

//Secant method for one function
//Parameters are function1, point 1, point 2, tolerance
csc350Lib_calculus_snle::SolutionNLE* Secant::solve(csc350Lib_calc_base::Function1D *f, float a, float b, float tol)
{
	//Solution to return
	csc350Lib_calculus_snle::SolutionNLE *nleSolution_ = new csc350Lib_calculus_snle::SolutionNLE;
	
	float fa = f->func(a);// function evaluated at a
	float fb = f->func(b);// function evaluated at b
	float d;
	int maxIterations = 100;
	int iterations=0;
	
	//Secant 
	for (int i = 2; i< maxIterations; i++) {
		//ensure that |fa| < |fb| to ensure convergence
		if (fabs(fa) > fabs(fb) ) {
			float temp = a;
			a = b;
			b = temp;
			fa = f->func(a);
			fb = f->func(b);
		}
		d = (b-a)/(fb-fa);
		b = a;
		fb = fa;
		d = d * fa;
		if (fabs(d) < tol ) {
			//upon convergence set and return solution
			cout << "convergance" << endl;
			nleSolution_->set_iterations(iterations);
			nleSolution_->set_func_value(f->func(a));
			nleSolution_->set_solution(a);
			return nleSolution_;
		}
		a = a - d;
		fa = f->func(a);
		iterations++;
	}
	//upon max iterations set and return solution
	nleSolution_->set_iterations(iterations);
	nleSolution_->set_func_value(f->func(a));
	nleSolution_->set_solution(a);
	return nleSolution_;
}//end solve function


//Secant Method to find where two functions intersect
//Parameters are function1, function2, point 1, point 2, tolerance
csc350Lib_calculus_snle::SolutionNLE* Secant::solve_for_2(csc350Lib_calc_base::Function1D *f,csc350Lib_calc_base::Function1D *f2, float a, float b, float tol)
{
	//solution to return
	csc350Lib_calculus_snle::SolutionNLE *nleSolution_ = new csc350Lib_calculus_snle::SolutionNLE;
	
	float fa = f->func(a) - f2->func(a);//function1 - function2 evaluated at a
	float fb = f->func(b) - f2->func(b);//function1 - function2 evaluated at b
	float d;
	int maxIterations = 20;
	int iterations;
	
	//Secant
	for (int i = 2; i< maxIterations; i++) {
		//check that |fa|< |fb| to ensure convergence
		if (fabs(fa) > fabs(fb) ) {
			float temp = a;
			a = b;
			b = temp;
			fa = f->func(a) - f2->func(a);
			fb = f->func(b) - f2->func(b);
		}
		d = (b-a)/(fb-fa);
		b = a;
		fb = fa;
		d = d * fa;
		if (fabs(d) < tol ) {
			//upon convergence set and return solution
			cout << "convergance" << endl;
			nleSolution_->set_iterations(iterations);
			nleSolution_->set_func_value(f->func(a)-f->func(b));
			nleSolution_->set_solution(a);
			return nleSolution_;
		}
		a = a - d;
		fa = f->func(a) - f2->func(a);
		iterations++;
	}
	//upon max iterations set and return solution
	nleSolution_->set_iterations(iterations);
	nleSolution_->set_func_value(f->func(a)-f->func(b));
	nleSolution_->set_solution(a);
	return nleSolution_;
}//end solve function

