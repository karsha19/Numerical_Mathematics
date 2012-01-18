#include <iostream>
#include "Newton.h"
using namespace std;

//Newton method for one function
//Parameters are function1, x point, max iterations, tolerance
csc350Lib_calculus_snle::SolutionNLE* Newton::solve(csc350Lib_calc_base::Function1D *f, float x, float max, float tol)
{
	//Solution to return
	csc350Lib_calculus_snle::SolutionNLE *nleSolution_ = new csc350Lib_calculus_snle::SolutionNLE;
	
	float  fx = f->func(x);     //f(x)
	float  dx;					//f'(x)
	int iterations = 1;
	float tolerance = 2 * tol;  //tolerance
	int maxIterations = 20;	//quit after max iterations
	
	for (int i =1; i<= maxIterations; i++) {
		dx = f->dfunc(x);//derivative of f at x
		
		//if the function converges set and return solution
		if (fabs(dx) < tolerance) {
			std::cout << "near root" << endl;
			nleSolution_->set_iterations(iterations);
			nleSolution_->set_func_value(f->func(x));
			nleSolution_->set_solution(x);
			return nleSolution_;
			
		}
		x = x - (fx / dx);
		fx = f->func(x);
	}
	//when max iterations are reached, set and return solution
	nleSolution_->set_iterations(iterations);
	nleSolution_->set_func_value(f->func(x));
	nleSolution_->set_solution(x);
	return nleSolution_;
	
}//end solve function


//Newton method to find where two functions intersect
//Parameters are function1, function2, x point, tolerance
csc350Lib_calculus_snle::SolutionNLE* Newton::solve_for_2(csc350Lib_calc_base::Function1D *f,csc350Lib_calc_base::Function1D *f2, float x, float tol)
{
	//Solution to return
	csc350Lib_calculus_snle::SolutionNLE *nleSolution_ = new csc350Lib_calculus_snle::SolutionNLE;
	
	float fx = f->func(x) - f2->func(x); // subtract one function from the other 
	float dx;
	int maxIterations = 10;
	int iterations = 1;
	float tolerance = 2* tol;
	
	
	for (int i =1; i<= maxIterations; i++) {
		dx = f->dfunc(x) - f2->dfunc(x);
		
		//if tolerance is reached set and return solution
		if (fabs(dx) < tolerance) {
			std::cout << "near root" << endl;
			nleSolution_->set_iterations(iterations);
			nleSolution_->set_func_value(f->func(x));
			nleSolution_->set_solution(x);
			return nleSolution_;
			
		}
		x = x - (fx / dx);
		fx = f->func(x) - f2->func(x);
		iterations++;
	}
	//when max iterations are reached, set and return solution
	nleSolution_->set_iterations(iterations);
	nleSolution_->set_func_value(f->func(x) - f2->func(x));
	nleSolution_->set_solution(x);
	return nleSolution_;
	
	
}//end solve_for_2
