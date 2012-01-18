#ifndef REALEIGENPAIR_H
#define REALEIGENPAIR_H
#include "Matrix.h"
using namespace csc350lib_linalg_base;

class RealEigenpair  {
public:
	RealEigenpair();
	RealEigenpair(float lamda, const Matrix<float> *vect);
	float get_eigenvalue(void) const;
	const Matrix<float> *get_eigenvector(void) const;
	void set(float, Matrix<float> *);
private:
	float eigenvalue_;
	const Matrix<float> *eigenvector_;
	
};



#endif