#include "RealEigenpair.h"

RealEigenpair::RealEigenpair(){};

//non-constant pointer to constant data meaning *vect cannot be modified
RealEigenpair::RealEigenpair(float lambda, const Matrix<float> *vect)
{
	eigenvalue_ = lambda;
	eigenvector_ = vect;
}


void RealEigenpair::set(float lambda, Matrix<float> *vect)
{
	eigenvalue_ = lambda;
	eigenvector_ = vect;
}
float RealEigenpair::get_eigenvalue(void) const
{
	return eigenvalue_;
}

const Matrix<float>* RealEigenpair::get_eigenvector(void) const
{
	return eigenvector_;
}