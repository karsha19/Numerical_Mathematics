#ifndef COLVEC_H
#define COLVEC_H
#include "Matrix.h"

namespace csc350lib_linalg_base {
	
template <class T> 
class ColVec : public Matrix<T>
{
	public:
		//Default constructor. Parameters: row = # of rows
		ColVec(int row): Matrix<T>(row, 1 ){};
		//Constructor. Parameters: row = # of rows, numbers = row entries
		ColVec(int row, T numbers[]) : Matrix<T>(row, 1, numbers) {};
};

	
}//end namespace
#endif