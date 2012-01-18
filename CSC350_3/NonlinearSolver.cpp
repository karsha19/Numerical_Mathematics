#include "NonlinearSolver.h"

namespace csc350Lib_calculus_snle{
	
	NonlinearSolver::NonlinearSolver(): nleSolution_(){}
	
	
	void NonlinearSolver::accept_nle_solution(SolutionNLE *solution)
	{
		nleSolution_ = solution;
	}
	
	
	float NonlinearSolver::get_solution_2()
	{	
		return nleSolution_->get_solution(); 
	}
	
	
	float NonlinearSolver::get_value_at_solution_2()
	{	
		return nleSolution_->get_value_at_solution();  
	}
	
	
	int NonlinearSolver::get_number_of_iterations_2()
	{	
		return nleSolution_->get_number_of_iterations();  
	}
	
	SolutionNLE* NonlinearSolver::solve_for_2(csc350Lib_calc_base::Function1D *f ,csc350Lib_calc_base::Function1D *f1, float x ,float y)
	{
		return	nleSolution_;
	}
	
	SolutionNLE* NonlinearSolver::solve_for_2(csc350Lib_calc_base::Function1D *f,csc350Lib_calc_base::Function1D *f2, float x, float, float tol)
	{
		return nleSolution_;
	}
	
	
	
	
}