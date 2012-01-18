#include "NonlinearSolver_bisection.h"

namespace csc350Lib_calculus_snle{
	
	NonlinearSolver_bisection::NonlinearSolver_bisection(void) : NonlinearSolver(){};
	
	//Solve a non linear equation using the bisection method
	//Retrun the answer as SolutionNle
	SolutionNLE* NonlinearSolver_bisection::solve(csc350Lib_calc_base::Function1D *f, float a, float b, float tol)
	{
		SolutionNLE *nleSolution_ = new SolutionNLE;
		float  fa = f->func(a);
		float  fb = f->func(b);
		float  c  = (a+b)/2;
		int iterations = 1;
		float tolerance = 2 * tol;
		while ((b-a > tolerance ) && (c != a) && (c != b)){
			//evalutae the function at the midpoint: c
			float fc = f->func(c);
			//if the function does not cross the y==0 line in the first half
			//then it must be the midpoint or in the second half
			//set the starting point to c(the old midpoint)
			if ((fc * fa) > 0 ) {
				a = c;
				fa = fc;
			}
			//check if c is the answer
			//if not, the answer must be in the first half
			else if (fc != 0){
				b = c;
				fb = fc;
			}
			else {
				nleSolution_->set_iterations(iterations);
				nleSolution_->set_func_value(f->func(c));
				nleSolution_->set_solution(c);
				return nleSolution_;
			}
		c = (a + b)/2;
		iterations++;


		}//end while loop
		 
		nleSolution_->set_iterations(iterations);
		nleSolution_->set_solution(f->func(c));
		nleSolution_->set_solution(c);
		return nleSolution_;
		
	 }//end solve function
	

}//end namespace