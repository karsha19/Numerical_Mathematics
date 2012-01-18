#ifndef ROWVEC_H
#define ROWVEC_H
#include "Matrix.h"

namespace csc350lib_linalg_base {
	
template <class T> 
class RowVec : public Matrix<T>
{
	public:
		//Default contructor. Parameters: col = # of columns
		RowVec(int col) : Matrix<T>(1, col){};
		//Constructor. Parameters: col = #of columns, numbers = column entries
		RowVec(int col, T numbers[] ) : Matrix<T>(1, col, numbers){};

};
	
}//end namespace


#endif