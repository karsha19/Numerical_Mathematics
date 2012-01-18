#include "EigensystemCalculator.h"

const RealEigenpair* EigensystemCalculator::get_largest_eigenpair( const Matrix<float> *m )
{
	Matrix<float>A(m->get_rows(), m->get_columns());
	A = *m;
	
	//max iterations should not exceed 100
	int max =100;
	
	float eigenValue =0;
	
	//boolean to determine sign of yMax * xMax
	bool is_sign_negative = false;
	
	//srand used to generate arbitrary matrix
	srand(time(0));
		
	//initialize arbitrary x vector
	float arbitrary_x_vector[A.get_rows()];
	for (int i =0; i<A.get_rows(); i++) {
		arbitrary_x_vector[i] = (rand()%10 ) +1;
	}
	
	//Eigenvector to be returned
	Matrix<float> *eigenVector;
	eigenVector = new Matrix<float>(A.get_rows(),1, arbitrary_x_vector);
	
	Matrix<float> y(A.get_rows(),1);
	
	int xIMax=0;
	int yIMax=0;
	float l = 0; //y.infinity_norm
	
	
	//normalized power iteration algorithm
	for (int i=0; i< max; i++) {
		
		y = A * (*eigenVector);
		
		l = y.infinity_norm();
		
		xIMax = (*eigenVector).get_max_i();
		
		yIMax = y.get_max_i();
		
		if ( ((*eigenVector).get_element(xIMax,0) * y.get_element(yIMax,0) ) <  0 ) {
			is_sign_negative = true;
		}
		else {
			is_sign_negative = false;
		}
		
		float oneover = 1/l;
		(*eigenVector) = y * oneover;
		
		
				
	}
	
	if (is_sign_negative) {
		eigenValue = -l;		
	}
	else {
		eigenValue = l;
	}
	

	float lambda =0;
		
	Matrix<float>xTranspose(1, (*eigenVector).get_rows() );

	Matrix<float>xtax(xTranspose.get_rows(), (*eigenVector).get_columns() );
	
	Matrix<float>xtx(xTranspose.get_rows(), (*eigenVector).get_columns() );
	
	Matrix<float>identity(A.get_rows(), A.get_columns());
	
	Matrix<float>ALambdaI(A.get_rows(), A.get_columns());
	
	Matrix<float>yy(A.get_rows(), 1);
	
	identity.set_identity();
	
	
	for (int i = 0; i < 10; i++) {
							
		xTranspose = (*eigenVector).transpose();
		
		xtax = xTranspose * ( A * (*eigenVector) );
		
		xtx = xTranspose * (*eigenVector);
		
		lambda = (xtax.get_element(0, 0)) / (xtx.get_element(0, 0));
		
		identity = identity * lambda;
		
		ALambdaI = A - identity;
		
		yy = ALambdaI.solve(*eigenVector);
		
		*eigenVector = yy	* (1/yy.infinity_norm());
		
	
	}
	
	const RealEigenpair *answer;
	answer = new RealEigenpair(eigenValue, eigenVector);
	return answer;
}

const RealEigenpair * EigensystemCalculator::get_eigenpairs(const Matrix<float> *m, int n)
{
	
	
	RealEigenpair * eigenPairs;
	eigenPairs = new RealEigenpair[n];
	
	Matrix<float>A(m->get_rows(), m->get_columns());
	A = *m;
	
	for (int k =0; k< n; k++) {
		
		//max iterations should not exceed 50
		int max =100;
		
		float eigenValue =0;
		
		//boolean to determine sign of yMax * xMax
		bool is_sign_negative = false;
		
		//srand used to generate arbitrary matrix
		srand(time(0));
		
		//initialize arbitrary x vector
		float arbitrary_x_vector[A.get_rows()];
		for (int i =0; i<A.get_rows(); i++) {
			arbitrary_x_vector[i] = (rand()%10 ) +1;
		}
		
		//Eigenvector to be returned
		Matrix<float> *eigenVector;
		eigenVector = new Matrix<float>(A.get_rows(),1, arbitrary_x_vector);
		
		Matrix<float> y(A.get_rows(),1);
		
		int xIMax=0;
		int yIMax=0;
		float l = 0;
		
		
		//normalized power iteration algorithm
		for (int i=0; i< max; i++) {
			
			y = A * (*eigenVector);
			
			l = y.infinity_norm();
			
			xIMax = (*eigenVector).get_max_i();
			
			yIMax = y.get_max_i();
			
			if ( ((*eigenVector).get_element(xIMax,0) * y.get_element(yIMax,0) ) <  0 ) {
				is_sign_negative = true;
			}
			else {
				is_sign_negative = false;
			}
			
			float oneover = 1/l;
			(*eigenVector) = y * oneover;
			
			
						
		}
		
		if (is_sign_negative) {
			eigenValue = -l;		
		}
		else {
			eigenValue = l;
		}
		
		
		float lambda =0;
		
		Matrix<float>xTranspose(1, (*eigenVector).get_rows() );
		
		Matrix<float>xtax(xTranspose.get_rows(), (*eigenVector).get_columns() );
		
		Matrix<float>xtx(xTranspose.get_rows(), (*eigenVector).get_columns() );
		
		Matrix<float>identity(A.get_rows(), A.get_columns());
		
		Matrix<float>ALambdaI(A.get_rows(), A.get_columns());
		
		Matrix<float>yy(A.get_rows(), 1);
		
		identity.set_identity();
		
		
		for (int i = 0; i < 10; i++) {
			
			xTranspose = (*eigenVector).transpose();
			
			xtax = xTranspose * ( A * (*eigenVector) );
			
			xtx = xTranspose * (*eigenVector);
			
			lambda = (xtax.get_element(0, 0)) / (xtx.get_element(0, 0));
			
			identity = identity * lambda;
			
			ALambdaI = A - identity;
			
			yy = ALambdaI.solve(*eigenVector);
			
			*eigenVector = yy	* (1/yy.infinity_norm());
			
			
		}
		
		eigenPairs[k].set(eigenValue, eigenVector);
		
		cout << "eigenValue= " << eigenValue << endl;
		cout << "eigenvector= " << endl;
		(*eigenVector).display();
		
		Matrix<float>u((*eigenVector).get_rows(),1);
		Matrix<float>uTranspose(1,(*eigenVector).get_rows());
		Matrix<float>ASubtractor(A.get_rows(), A.get_columns());
		
		u = (*eigenVector) * (eigenValue * (  1 / (pow( (*eigenVector).euclidean_norm(), 2)  )));
		
		uTranspose = u.transpose();
		
		ASubtractor = (*eigenVector) * uTranspose;
		
		A = A - ASubtractor;
		
		
	}//end for loop
	
	return eigenPairs;
	
}// end 2
