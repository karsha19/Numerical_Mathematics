#include <iostream>
#include "Matrix.h"
#include "NonlinearSolver_bisection.h"
#include "Function1D.h"
#include "EigensystemCalculator.h"
#include "RealEigenpair.h"
#include "Assignment3.h"
#include <stdlib.h>

using namespace std;
using namespace csc350Lib_calc_base;
using namespace csc350Lib_calculus_snle;
using namespace csc350lib_linalg_base;



int main(int argc, char** argv) {
	
	struct point {
		float x,
			  y,
			  z, 
			  distance;
	};
	
	//random points I made up 
	point bowl[]={	 
					 0  , 0  ,  1,  0,
					 0  ,-1  ,  1,  0,
					 .5 ,-.5 ,  1,  0,
					 1  , 0  ,  1,  0,
					 .5 , .5 ,  1,  0,
					  0 , 1  ,  1,  0,
					-.5 ,.5  ,  1,  0,
					-1  , 0  ,  1,  0,
					-.5 ,-.5 ,  1,  0,
					 0  , -2 ,  2,  0,
					2.5 ,-2.5,  2,  0,
					2   , 0  ,  2,  0,
					2.5 , 2.5,  2,  0,
					0   , 2  ,  2,  0,
					-2.5, 2.5,  2,  0,
					-2  , 0  ,  2,  0,
					-2.5,-2.5,  2,  0,
					0   ,-3  , 2.5, 0,
					3.5 ,-3.5, 2.5, 0,
					3   , 0  , 2.5, 0,
					3.5 , 3.5, 2.5, 0,
					0   , 3  , 2.5, 0,
					-3.5, 3.5, 2.5, 0,
					-3  , 0  , 2.5, 0,
					-3.5,-3.5, 2.5, 0,
					0   , -4 , 3,   0,
					4.5 ,-4.5, 3,   0,
					4   , 0  , 3,   0,
					4.5 , 4.5, 3,   0,
					0   , 4  , 3,   0,
					-4.5, 4.5, 3,   0,
					-4  , 0  , 3,   0,
					-4.5,-4.5, 3,   0,
					0   ,-5  , 3.5, 0,
					5.5 ,-5.5, 3.5, 0,
					5   , 0  , 3.5, 0,
					5.5 , 5.5, 3.5, 0,
					0   , 5  , 3.5, 0,
					-5.5,5.5 , 3.5, 0,
					-5  , 0  , 3.5, 0,
					-5.5,-5.5, 3.5, 0 };
	
	//choice is the point i choose to find the curvature(eigenvalues) of
	point choice;
	choice.x =-2;
	choice.y =0;
	choice.z =2;
	
	//numberofpoints = 15 points closest to choice
	int numberOfPoints = 15;
	for (int i = 0; i <sizeof(bowl)/sizeof(point); i++) {
		
		bowl[i].distance = sqrt(pow(bowl[i].x -choice.x , 2) +pow(bowl[i].y -choice.y, 2) +pow(bowl[i].z -choice.z, 2));
	}
	
	for (int i=numberOfPoints; i < sizeof(bowl)/sizeof(point); i++) {
		
		for (int j=0; j< numberOfPoints; j++) {
			if (bowl[i].distance < bowl[j].distance) {
				
				point temp;
				temp.x = bowl[j].x;
				temp.y = bowl[j].y;
				temp.z = bowl[j].z;
				temp.distance = bowl[j].distance;
				
				bowl[j].x = bowl[i].x;
				bowl[j].y = bowl[i].y;
				bowl[j].z = bowl[i].z;
				bowl[j].distance = bowl[i].distance;
				
				bowl[i].x = temp.x;
				bowl[i].y = temp.y;
				bowl[i].z = temp.z;
				bowl[i].distance = temp.distance;
			}
			
			
		}
	}
	//15x6 matrix of points closest to choice
	Matrix<float> points(numberOfPoints,6);
	
	//15x1 column vector of z points
	Matrix<float> zColumn(numberOfPoints,1);
	
	//set up matricies to solve equation 
	// x^2 + 2xy + y^2 + x + y + 1  = z
	for (int i=0; i< numberOfPoints; i++) {
		points.set_element(i, 0, pow(bowl[i].x, 2));
		points.set_element(i, 1, (2*bowl[i].x*bowl[i].y));
		points.set_element(i, 2, pow(bowl[i].y, 2));
		points.set_element(i, 3, bowl[i].x);
		points.set_element(i, 4, bowl[i].y);
		points.set_element(i, 5, 1);
		
		zColumn.set_element(i,0, bowl[i].z);
		
	}
	
	Matrix<float> answer(6,1);
	answer = points.solve_w_householder(zColumn);
	
	
	Assignment3 curve;
	curve.set_coefficients(answer.get_element(0, 0),answer.get_element(1, 0),answer.get_element(2, 0),
						   answer.get_element(3, 0),answer.get_element(4, 0),answer.get_element(5, 0));
	
	float curvature = curve.compute_curvature(choice.x, choice.y);
	
	cout << "curvature is " << curvature << endl;
	
	//testing & debugging	
	/*float eigenarray[9] = {-4,-4,-10,
						   -6,0,3,
						   1,1,3};
	
	
	const Matrix<float> *eigenmatrix;
	eigenmatrix = new Matrix<float>(3,3, eigenarray);
	
	EigensystemCalculator *Ecalculator;
	Ecalculator = new EigensystemCalculator();
	
	const RealEigenpair *eigenpair;
	eigenpair= new RealEigenpair();
	eigenpair = Ecalculator->getLargestEigenpair(eigenmatrix);
	
	const RealEigenpair * pairs;
	pairs = new RealEigenpair[3];
	
	pairs = Ecalculator->getEigenpairs(eigenmatrix,3);
	
	cout << "The eigenvalue is " << eigenpair->getEigenvalue() << endl;
	cout << "The eigenvector is "<< endl;
	eigenpair->getEigenvector()->display();
	
	float hholderarray[24] = {	1,1,1,1,
								1,0,1,1,
								0,1,1,1,
								0,1,0,1,
								1,1,1,1,
								0,0,1,1 };
	float zzColumn[6] ={16,10,11,8,16,5};
	
	Matrix<float> hholder(6,4, hholderarray);
	Matrix<float> b(6,1,zzColumn);
	
	

	hholder.solve_w_householder(b);*/
		
	return 0;

}


