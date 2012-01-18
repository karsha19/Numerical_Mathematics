#ifndef NONLINEARSOLVER_BISECTION_H
#define NONLINEARSOLVER_BISECTION_H
#include "NonlinearSolver.h"

namespace csc350Lib_calculus_snle {

	//Solve a non linear equation using the bisection method
	//This subclass simply implements the solve method of the parent class.
	class NonlinearSolver_bisection : public NonlinearSolver {
		public:
			NonlinearSolver_bisection(void);                    
			~NonlinearSolver_bisection(void){};  
			SolutionNLE* solve(csc350Lib_calc_base::Function1D *f, float a, float b, float tol);
	};
	
}
#endif