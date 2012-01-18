#ifndef NONLINEARSOLVER_H
#define NONLINEARSOLVER_H
#include "SolutionNLE.h"
#include "NonlinearSolver.h"

namespace csc350Lib_calculus_snle{
	/*
	 The pure abstract class NonlinearSolver
	 This is the parent class for the various solvers that you will implement. 
	 This class and its subclasses will operate on Function1D objects.
	*/
	class NonlinearSolver {
	public:
		//Default constructor
		NonlinearSolver(void);
		
		//Destructor
		virtual ~NonlinearSolver(void){};
		
		/*This is a pure abstract method that should be implemented by subclasses. 
		 The parameters are a pointer to the NLE's function, the endpoints of the search's bracket, 
		 and the tolerance of the search. The function (again, implemented by the subclasses) 
		 should return a pointer to a SolutionNLE object.
		 */
		virtual SolutionNLE* solve(csc350Lib_calc_base::Function1D *f, float a, float b, float tol) =0;
		
		//used by the Newton method solver child class
		virtual SolutionNLE* solve_for_2(csc350Lib_calc_base::Function1D *f,csc350Lib_calc_base::Function1D *f2, float, float);
		
		//used by the Secant method solver child class
		virtual SolutionNLE* solve_for_2(csc350Lib_calc_base::Function1D *f,csc350Lib_calc_base::Function1D *f2, float x, float, float tol);
		
		virtual void  accept_nle_solution(SolutionNLE * solution);
		
		//Getter methods
		virtual float get_solution_2();
		virtual float get_value_at_solution_2();
		virtual int   get_number_of_iterations_2();
	private:
		SolutionNLE *nleSolution_;
		
		
	};//end nonlinearsolver class
	
	
	
}//end calculus_snle namespace

#endif