#ifndef SOLUTIONNLE_H
#define SOLUTIONNLE_H
#include "Function1D.h"

namespace csc350Lib_calculus_snle {
	/*
	 The result of an NLE search problem is returned as an instance of the SolutionNLE class.
	 An element of the class SolutionNLE stores as its instance variables
	 the value of solution estimate x when the search stopped,
	 the value of f(x) when the search stopped,
	 the number of iterations when the search stopped.
	*/
	class SolutionNLE{
		
	public:
		//Default constructor
		SolutionNLE(){}; 
		
		//Return the value of x when the search is stopped
		float get_solution(void);
		
		//Return the value of f(x) when the search is stopped
		float get_value_at_solution(void);
		
		//Return the number of iterations when the search is stopped
		int   get_number_of_iterations(void) const; 
		
		//Setter methods to set class members
		void  set_solution(float);
		void  set_func_value(float);
		void  set_iterations(int);
		~SolutionNLE();
	private:
		
		//The value of x when the search is stopped
		float solution_;
		
		//The value of f(x) when the search is stopped
		float funcValue_;
		
		//The number of iterations when the search is stopped
		int   numOfIterations_;
		
	};//end SoultionNLE class
	
}//end calculus_snle namespace

#endif