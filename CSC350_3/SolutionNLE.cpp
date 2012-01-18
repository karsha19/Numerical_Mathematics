#include "SolutionNLE.h"

namespace csc350Lib_calculus_snle {
	
	float SolutionNLE::get_solution(void)
	{
		return solution_;
	}
	
	float SolutionNLE::get_value_at_solution(void)
	{
		return funcValue_;
	}
	
	int SolutionNLE::get_number_of_iterations(void) const
	{
		return numOfIterations_;
	}
	
	void SolutionNLE::set_solution(float answer)
	{
		solution_ = answer;
	}
	
	void SolutionNLE::set_iterations(int iterations)
	{
		numOfIterations_ = iterations;
	}
	
	void SolutionNLE::set_func_value(float answer)
	{
		funcValue_ = answer;
	}
	
	SolutionNLE::~SolutionNLE(){}
	
}