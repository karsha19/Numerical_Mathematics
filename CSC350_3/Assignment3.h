#ifndef ASSIGNMENT3_H
#define ASSIGNMENT3_H
#include "Matrix.h"
#include "EigensystemCalculator.h"
#include "RealEigenpair.h"

class Assignment3 {

public:
	Assignment3();
	~Assignment3(){};
	float dot_j(int);
	void set_coefficients(float, float, float, float, float, float);
	void compute_hessian();
	void compute_gradient(float,float);
	void compute_gradient_norm(float, float);
	void compute_c_matrix();
	void compute_c_eigenvalues(float *);
	float compute_curvature(float, float);
	
private:
	
	float	a_,b_,c_,d_,e_,f_,
			gradientNormal_,
			gradientNSquared_;

	float	gradient_[3];
	Matrix<float> hessian_;
	Matrix<float> cMatrix_;
	
};

Assignment3::Assignment3(): hessian_(3,3), cMatrix_(3,3){};
//compute dot j
float Assignment3::dot_j(int j)
{
	float item;
	item = gradient_[0] * hessian_.get_element(j,0) + gradient_[1] * hessian_.get_element(j,1);
	return item;
}

void  Assignment3::set_coefficients( float a, float b, float c, float d, float e, float f)
{
	a_ = a;
	b_ = b;
	c_ = c;
	d_ = d;
	e_ = e;
	f_ = f;
	
}

void Assignment3::compute_hessian()
{
	//compute hessian matrix
	hessian_.set_element(0, 0, 2*a_);
	hessian_.set_element(0, 1, 2*b_);
	hessian_.set_element(1, 0, 2*b_);
	hessian_.set_element(1, 1, 2*c_);
}
void Assignment3::compute_gradient(float choiceX, float choiceY)
{
	gradient_[0] = 2*a_*choiceX+2*b_*choiceY+d_;
	gradient_[1] = 2*b_*choiceX+2*c_*choiceY+e_;
	gradient_[2] = 1;
	
}
void Assignment3::compute_gradient_norm(float choiceX, float choiceY)
{
	//compute gradientNormal
	float gradientNorm_[3];
	gradientNorm_[0] = pow(2*a_*choiceX+2*b_*choiceY+d_,2);
	gradientNorm_[1] = pow(2*b_*choiceX+2*c_*choiceY+e_,2);
	gradientNorm_[2] = 1;
	gradientNormal_ = sqrt(gradientNorm_[0] + gradientNorm_[1] + gradientNorm_[2]);
	gradientNSquared_ = pow(gradientNormal_,2);
}

void Assignment3::compute_c_matrix()
{
	//compute c matrix
	for (int i=0; i < 3; i++) {
		for (int j=0; j<3; j++) {
			cMatrix_.set_element(i, j, ((hessian_.get_element(i,j)*gradientNormal_) - ((gradient_[i] * dot_j(j))/gradientNormal_))/gradientNSquared_);
		}
	}
	cout << "gradient 1"<< endl;
		cout << gradient_[0]<< endl;
	cout << "gradient 2"<< endl;
		cout << gradient_[1]<< endl;
	cout << "gradient 3"<< endl;
		cout << gradient_[2]<< endl<< endl;
	
	cout << "hessian " << endl;
	hessian_.display();
	
	cout << "cMatrix " << endl;
	cMatrix_.display();
}

void Assignment3::compute_c_eigenvalues(float *k)
{
	EigensystemCalculator *eigensystemCalc;
	eigensystemCalc = new EigensystemCalculator();
	
	const RealEigenpair * pairs;
	pairs = new RealEigenpair[2];
	
	pairs = eigensystemCalc->get_eigenpairs(&cMatrix_,2);
	
	k[0] = pairs[0].get_eigenvalue();
	k[1] = pairs[1].get_eigenvalue();
}

float Assignment3::compute_curvature(float x, float y)
{
	compute_hessian();
	
	compute_gradient(x, y);
	
	compute_gradient_norm(x, y);
	
	compute_c_matrix();
	
	float k[2];
	
	compute_c_eigenvalues(k);
	
	float curvature = k[0] * k[1];
	
	return curvature;
}

#endif


