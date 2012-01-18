#define PI 3.14159
#include "Interpolation.h"

using namespace csc350lib_linalg_base;

Interpolation::Interpolation(){}

PolyFunction1D Interpolation::polynomial_case(int numberOfPoints, float inputPoints[][2] )
{
			
	//number of matrix entries
	int number_of_entries = numberOfPoints * numberOfPoints;
	
	//array to initialize b column
	float tempArray[numberOfPoints];
	for (int i = 0; i < numberOfPoints; i++) {
		tempArray[i] = inputPoints[i][1];
	}
	
	//setup bcolumn as y coordinates
	ColVec<float> bColumn(numberOfPoints,tempArray);
	
	// define array to initialize matrix
	float matrixNumbers[number_of_entries];
	
	//initialize matrix_array
	int count = 0;
	for (int i = 0; i< numberOfPoints; i++) {
		for (int j = 0; j< numberOfPoints; j++) {
			matrixNumbers[count++] = pow(inputPoints[i][0],j);
		}
	}
	
	//initialize actual matrix using matrix_array
	Matrix<float> temp(numberOfPoints,numberOfPoints, matrixNumbers);
	
	//define matrix to hold answer
	Matrix<float> answer(numberOfPoints,1);
	answer = temp.solve( bColumn);
	
	//define array to hold polyomial coefficients
	float coefficients[numberOfPoints];
	for (int i =0; i< numberOfPoints; i++) {
		coefficients[i] = answer.get_element(i, 0);
	}
	
	PolyFunction1D tempf(numberOfPoints,coefficients);
	return tempf;
	


}


float Interpolation::lagrange_case(int numberOfPoints, float inputPoints[][2], float x )
{
	//number of points
	int n = numberOfPoints;
	
	float fx = 0;//the value of f(x)
	
	for(int i=0; i<n; i++){
		float lagrange = 1;
		for(int j=0; j<n; j++){
			if(i != j){
				lagrange =lagrange * (x - inputPoints[j][0])/(inputPoints[i][0]-inputPoints[j][0]);
				
			}
		}
		fx = fx + lagrange * inputPoints[i][1];
		
	}
	
	return fx;
	
	
}

PolyFunction1D Interpolation::trigonometric_case( int numberOfPoints, float inputPoints[][2], float L, float n )
{
	
	int dimension =  2 * n + 1;
	
	int total_entries = pow(dimension, 2.0) ;
	
	float matrix_array[total_entries];
	float bcolumn_array[dimension];
	
	
	float x;
	
	int count = 0;
	
	for (int rows = 0; rows < dimension; rows++) {
		
		x = inputPoints[rows][0];
		
		bcolumn_array[rows] = inputPoints[rows][1];// y coordinates
		
		
		matrix_array[count++] = inputPoints[rows][1]/dimension;// k in the function?
		
		for (int columns = 1; columns <= n; columns++) {
			matrix_array[count++] = cos((2 * columns * PI * x)/ L );
			
		}
		for (int columns = 1; columns <= n; columns++) {
			matrix_array[count++] = sin((2 * columns * PI * x)/ L );
					}
		
	}
	
	
	Matrix<float> temp(dimension,dimension, matrix_array);
	temp.display();
	Matrix<float> bcolumn(dimension, 1, bcolumn_array);
	Matrix<float> coefficients_matrix(dimension,1);
	coefficients_matrix = temp.solve( bcolumn);
 	coefficients_matrix.display();
	
	
	//define array to hold polyomial coefficients
	float coefficients[numberOfPoints];
	
	//assign coefficients
	for (int i =0; i< numberOfPoints; i++) {
		coefficients[i] = coefficients_matrix.get_element(i, 0);
	}
	
	PolyFunction1D tempf(dimension, coefficients);
	
	return tempf;
	
}

