#ifndef MATRIX_H
#define MATRIX_H
#include <vector>
#include <iostream>
#include <cmath>
#include <assert.h>
using namespace std;

namespace csc350lib_linalg_base {
	
template <typename T>
class Matrix{	
public:
	//constructors
	Matrix(int row, int col, T numbers[] );
	Matrix(int row, int col);
	//overload operators to work on matricies
	const Matrix<T>& operator=(const Matrix<T> & ); 
	const Matrix<T>  operator+(const Matrix<T> & ); 
	const Matrix<T>  operator-(const Matrix<T> & );
	const Matrix<T>  operator*(float s);		  
	const Matrix<T>  operator*(const Matrix<T> & ); 
	//transpose a matrix
	const Matrix<T>  transpose();
	//apply back substitution to a triangluar matrix
	const Matrix<T>  back_substitution( const Matrix<T> &);
	const Matrix<T>  forward_substitution( const Matrix<T> &, int []);
	//factorize a matrix using the crout method without a pivot 
	const Matrix<T>  crout_no_pivot();
	//factorize a matrix using the crout method with a pivot
	const Matrix<T>  crout_pivot(int []);
	//factorize a matrix using the householder method
	void householder();
	//set a matrix to the identity matrix
	const Matrix<T>  set_identity( );
	//solve a matrix system of Ax=b.  
	//A is the calling matrix b is the parameter matrix
	const Matrix<T>  solve(Matrix<T> & );
	
	const Matrix<T> solve_w_householder(Matrix<T> & );
	
	T     get_element(int i, int j);
	T     get_matrix_determinant();
	T     infinity_norm();
	T     euclidean_norm();
	T     get_max_i();
	int   get_rows()const;
	int   get_columns()const;
	void  set_element(int i, int j, T );
	void  display()const;
	void  set_rows( int );
	void  set_columns( int );
	void  replace_w_zeros();
	T  get_matrix_entry(int , int );
	void  init_vector(T numbers[]);
	
	
	
	
private:
	//vector that will act as matrix
	vector<vector<T> > objMatrix;
	
	// rows of matrix
	int rows_; 
	
	// columns of matrix
	int columns_;
	
};//-----------------------------------
	
	//constructor without initialzation
	template <typename T>
	Matrix<T>::Matrix(int row, int col)
	{
		rows_ = row;   //set rows_
		columns_ = col;//set columns_
		
		// configure vector into matrix arrangement
		objMatrix.resize(row);
		for (int i = 0; i < row; ++i)
			objMatrix[i].resize(col);
	}


	//constructor with initializatoin
	template <typename T>
	Matrix<T>::Matrix(int row, int col, T numbers[] )
	{
		//set rows_
		rows_ = row; 
		
		//set columns_
		columns_ = col;
		
		// configure vector into matrix arrangement
		objMatrix.resize(row);
		for (int i = 0; i < row; ++i)
			objMatrix[i].resize(col);
		
		//input given numbers into matrix
		init_vector( numbers );
	}

	//init_vector called from constructor w/ initilization 
	template <typename T>
	void Matrix<T>::init_vector( T numbers[] )
	{
		int count = 0;
		for (int i = 0; i < rows_; i++)
			for (int j =0; j< columns_; j++)
				objMatrix[i][j] = numbers[count++];
	}
	
	//Return matrix element at [i][j]
	template <typename T>
	T Matrix<T>::get_element(int i, int j)
	{	return objMatrix[i][j]; }
	
	//set matrix element at [i][j]
	template <typename T>
	void Matrix<T>::set_element(int i, int j, T val)
	{	objMatrix[i][j] = val; }
	
	//set number of rows_ of matrix
	template <typename T>
	void Matrix<T>::set_rows(int r)
	{ rows_ =r; }

	//set number of columns_ of matrix
	template <typename T>
	void Matrix<T>::set_columns(int c)
	{ columns_ = c;}

	//return rows_ of matrix
	template <typename T>
	int Matrix<T>::get_rows()const
	{ return rows_;}

	//return columns_ of matrix
	template <typename T>	
	int Matrix<T>::get_columns()const
	{ return columns_;}

	//return desired matrix entry
	template <typename T>
	T Matrix<T>::get_matrix_entry(int a, int b)
	{	return objMatrix[a][b]; }
	
	
	//find determinant of matrix
	template <typename T>
	T Matrix<T>::get_matrix_determinant()
	{
		Matrix<T> temp( (*this).rows_, (*this).columns_ );
		temp = (*this).crout_no_pivot();
		T determinant = 1;
		for (int i = 0; i< (*this).rows_; i++) {
			determinant = determinant * (*this).objMatrix[i][i];
		}
		return determinant;
	}

	
	//assignment operator
	template <typename T>
	const Matrix<T>& Matrix<T>::operator = (const Matrix<T>& m){
		//assert(m.rows_ != rows_ && m.columns_ != columns_ );
	
		for(unsigned i=0; i<m.rows_; i++)
			for(unsigned j=0; j<m.columns_; j++)
				objMatrix[i][j] = m.objMatrix[i][j];
		
		return *this;
	}

	// + operator
	template <typename T>
	const Matrix<T> Matrix<T>::operator+( const Matrix<T>& m) { 
		//assert((*this).columns_==m.columns_ && (*this).rows_==m.rows_); 
		Matrix<T> temp(rows_,columns_); 
		for (int r=0;r<rows_;r++) 
			for (int c=0;c<columns_;c++) 
				temp.objMatrix[r][c]=objMatrix[r][c]+m.objMatrix[r][c]; 
		return temp; 
	} 
	// - operator
	template <typename T>
	const Matrix<T> Matrix<T>::operator-( const Matrix<T>& m) { 
		//assert((*this).columns_==m.columns_ && (*this).rows_==m.rows_); 
		Matrix<T> temp(rows_,columns_); 
		for (int r=0;r<rows_;r++) 
			for (int c=0;c<columns_;c++) 
				temp.objMatrix[r][c]=objMatrix[r][c]-m.objMatrix[r][c]; 
		return temp; 
	} 

	//multiply matrix by a scalar
	template <typename T>
	const Matrix<T> Matrix<T>::operator*( float s) { 
		Matrix<T> temp(rows_,columns_); 
		for (int r=0;r<rows_;r++) 
			for (int c=0;c<columns_;c++) 
				temp.objMatrix[r][c]=objMatrix[r][c]*s; 
		return temp; 
	} 

	// * operator
	template <typename T>
	const Matrix<T> Matrix<T>::operator*(const Matrix<T>& m){ 
		//assert((*this).columns_==m.rows_); 
		Matrix<T> temp(rows_,m.columns_); 
		for (int r=0;r<rows_;r++) { 
			for (int c=0;c<m.columns_;c++) { 
				for (int i=0;i<columns_;i++) { 
					temp.objMatrix[r][c]+=objMatrix[r][i]*m.objMatrix[i][c]; 
				} 
			} 
		} 
		return temp; 
	}


	// transpose matrix
	template <typename T>
	const Matrix<T> Matrix<T>::transpose()
	{
		Matrix<T> temp((*this).columns_,(*this).rows_);
		
		for (int i=0; i < (*this).rows_; i++)
			for (int j=0; j < (*this).columns_; j++)
			{
				//T x = (*this).objMatrix[i][j];
				temp.objMatrix[j][i] = (*this).objMatrix[i][j];
			}
		return temp;
	}

	//output matrix to console
	template <typename T>
	void Matrix<T>::display()const
	{
		for(int i=0; i<(*this).rows_; i++)
			for(int j=0; j<(*this).columns_; j++)
			{
				if (j == (*this).columns_ -1){
					cout << (*this).objMatrix[i][j] << endl;}
				else {
					cout << (*this).objMatrix[i][j] << "," ;}
					
			}
		cout << endl;
	}
	
	
	template <typename T>
	const Matrix<T> Matrix<T>::set_identity()
	{
		for (int i =0; i< (*this).rows_; i++)
			(*this).objMatrix[i][i]= 1;
		
		return *this;
	}
		


	//backsubstitution
	template <typename T>
	const Matrix<T> Matrix<T>::back_substitution( const Matrix<T>& m2)
	{
		Matrix<T> temp((*this).rows_, 1);
		int n = (*this).rows_-1;
		float sum;

		temp.objMatrix[n][0] = m2.objMatrix[n][0] / (*this).objMatrix[n][n] ;
		for (int i = n - 1; i >=0; i--) {
			sum = m2.objMatrix[i][0];
			for (int j = i+1; j <= n; j++) {
				sum = sum - ((*this).objMatrix[i][j] * temp.objMatrix[j][0]);
				
			}
			
			temp.objMatrix[i][0] = sum/(*this).objMatrix[i][i];
			
		}
		return temp;
	}


	//forward substitution
	template <typename T>
	const Matrix<T> Matrix<T>::forward_substitution( const Matrix<T>& m2, int pArray[])
	{
		
		Matrix<T> temp((*this).rows_, 1);
		int n = (*this).rows_-1;
		float sum;
		
		temp.objMatrix[pArray[0]][0] = m2.objMatrix[pArray[0]][0] / 1.0;
		for (int i = 1; i <=n; i++) {
			sum = m2.objMatrix[pArray[i]][0];
			for (int j=0; j <= i-1; j++) {
				sum = sum - ((*this).objMatrix[pArray[i]][j] * temp.objMatrix[pArray[j]][0]);
			}
			
			temp.objMatrix[pArray[i]][0] = sum/ 1.0;
		}
		return temp;
		
	}

	// crout LU factorization without pivot
	template <typename T>
	const Matrix<T> Matrix<T>::crout_no_pivot()
	{
		Matrix<T> temp((*this).rows_,(*this).rows_);
		Matrix<T> lower((*this).rows_,(*this).rows_);
		Matrix<T> upper((*this).rows_,(*this).rows_);
		Matrix<T> product((*this).rows_,(*this).rows_);
		int n = (*this).rows_;
		float sum;
		for (int j = 0; j<n; j++){
			for (int i = 0; i<=j; i++){
				sum = (*this).objMatrix[i][j];
				for (int k=0; k <= i; k++) {
					sum = sum - (temp.objMatrix[i][k] * temp.objMatrix[k][j]);
				}
				temp.objMatrix[i][j]= sum;
				upper.objMatrix[i][j]=sum;
			}
			for (int i =j+1; i<n; i++){
				sum = (*this).objMatrix[i][j];
				for (int k=0; k <= j; k++) {
					sum = sum - (temp.objMatrix[i][k] * temp.objMatrix[k][j]);
				}
				temp.objMatrix[i][j]= sum/temp.objMatrix[j][j];
				lower.objMatrix[i][j] = sum/temp.objMatrix[j][j];
			}
		}
		
		return temp;						 
				
	}

	
	template <typename T>
	const Matrix<T> Matrix<T>::crout_pivot(int pArray[] )
	{
		Matrix<T> upper((*this).rows_,(*this).rows_);
		Matrix<T> lower((*this).rows_,(*this).rows_);
		Matrix<T> multi((*this).rows_,(*this).rows_);
		Matrix<T> temp((*this).rows_,(*this).rows_);
		int n = (*this).rows_ -1;

		int ratioindex=0;
		float ratio = 0;
		float temporary;//debugging
		
		//find largest column in each row
		float largest_in_row[(*this).rows_];
		for (int i =0; i < (*this).rows_; i++) {
			largest_in_row[i] = 0;
			for (int j =0; j < (*this).columns_; j++) {
				if (abs((*this).objMatrix[i][j]) > largest_in_row[i]) {
					largest_in_row[i] = abs((*this).objMatrix[i][j]);
					temporary = largest_in_row[i]; //debugging
				}
			}
		}
		
		
		
		float sum;
		for (int j = 0; j<=n; j++){
			//upper
			for (int i = 0; i<=j-1; i++){
				sum = (*this).objMatrix[pArray[i]][j];
				for (int k=0; k <=i-1; k++) {
					sum = sum - (lower.objMatrix[pArray[i]][k] * upper.objMatrix[pArray[k]][j]);
				}
				upper.objMatrix[pArray[i]][j]= sum;
				temp.objMatrix[pArray[i]][j]= sum;
				temporary = sum;
				upper.display();
			}
			//multi
			for (int i =j; i<=n; i++){
				sum = (*this).objMatrix[pArray[i]][j];
				for (int k=0; k <=j-1; k++) {
					sum = sum - (lower.objMatrix[pArray[i]][k] * upper.objMatrix[pArray[k]][j]);
				}
				multi.objMatrix[pArray[i]][j]= sum;
				multi.display();
				temporary = sum;
			}
			
			temporary = multi.objMatrix[pArray[j]][j];
			temporary = largest_in_row[pArray[j]];
			
			ratio = 0;
			temporary = ratio;
			for (int i = j; i<=n; i++) {
				temporary = multi.objMatrix[pArray[i]][j];
				temporary = largest_in_row[pArray[i]];
				if (abs(multi.objMatrix[pArray[i]][j])/largest_in_row[pArray[i]] > ratio) {
					ratio = abs(multi.objMatrix[pArray[i]][j])/largest_in_row[pArray[i]];
					temporary = ratio;
					ratioindex = pArray[i];
				}
			}
			//swap if not largest
			if (ratioindex != j) {
				cout << "we in here" << endl;
				cout << pArray[0] << endl;
				cout << pArray[1] << endl;
				swap(pArray[ratioindex], pArray[j]);
				cout << pArray[0] << endl;
				cout << pArray[1] << endl;
				swap(largest_in_row[pArray[ratioindex]], largest_in_row[pArray[j]]);
				temporary = pArray[ratioindex];
			}
			temp.objMatrix[pArray[j]][j] = multi.objMatrix[pArray[ratioindex]][j];
			upper.objMatrix[pArray[j]][j] = multi.objMatrix[pArray[ratioindex]][j];
			upper.display();
			//multiplier
			float multiplier = 0 ;
			if (upper.objMatrix[pArray[j]][j] != 0) {
				temporary = upper.objMatrix[pArray[j]][j];
				multiplier = 1/upper.objMatrix[pArray[j]][j];
			}
			//lowermatrix
			cout << multiplier << " multiplier for j = " << j << endl;
			for (int i = j+1; i <=n; i++) {
				temporary = multi.objMatrix[pArray[i]][j];
				temporary = multiplier;
				lower.objMatrix[pArray[i]][j] = multiplier * multi.objMatrix[pArray[i]][j];
				temp.objMatrix[pArray[i]][j]= multiplier * multi.objMatrix[pArray[i]][j];
				temporary = multiplier * multi.objMatrix[pArray[i]][j];
			}
			
			
		}
		
		for (int i = 0 ; i<(*this).rows_; i++) {
			cout << pArray[i] << "parray"<< endl;
			cout << largest_in_row[i] << "largest_in_row"<<endl;
		}
		upper.display();
		lower.display();
		multi.display();
		temp.display();
		
		return temp;						 
		
	}

	// SLE solve method
	template <typename T>
	const Matrix<T> Matrix<T>::solve( Matrix<T>& m2)
	{
		//init permutation array
		int numbers[(*this).rows_];
		for (int i =0; i < (*this).rows_; i++) {
			numbers[i] = i;
		}
		
		Matrix<T> temp((*this).rows_,1);
		Matrix<T> answer((*this).rows_,1);
		Matrix<T> factorizedMatrix((*this).rows_,(*this).columns_);
		factorizedMatrix = (*this).crout_no_pivot();
		
		temp = factorizedMatrix.forward_substitution(m2, numbers);
		
		answer = factorizedMatrix.back_substitution(temp);
		
		return answer;
		
		
	}
	
	template <typename T>
	const Matrix<T> Matrix<T>::solve_w_householder( Matrix<T> & bColVector)
	{
		
		int inputRow = (*this).get_rows();
		int inputColumn = (*this).get_columns();
		(*this).display();
		Matrix<T> qRMatrix(inputRow,inputRow);
		qRMatrix.set_identity();
		
		Matrix<T> h(inputRow, inputRow);
		for(int j=0; j<inputColumn  ; j++){
			//a e v h, same vectors/matricies from algorithm in notes
			Matrix<T> a(inputRow, 1);
			Matrix<T> e(inputRow, 1);
			Matrix<T> v(inputRow, 1);
			
			
			//create identity matrix 
			Matrix<T> identity(inputRow,inputRow);
			identity.set_identity();
			
			
			int rowIndex = j;
			
			//make a = to current column
			for(int i=j; i<inputRow; i++){
				a.objMatrix[i][0] = (*this).objMatrix[rowIndex++][j]; 
			}
				//cout<< "a= " << endl;
				//a.display();
			
			//initialize the e colomn vector
			e.objMatrix[j][0]=1;
			
				//cout << "e= " << endl;
				//e.display();
			
			//compute norm of current column
			float elementSqSum =0;
			for(int i=j; i<inputRow; i++){
				elementSqSum += pow(a.objMatrix[i][0],2);
			}
			float aNorm = sqrt(elementSqSum);
				//cout << "anorm= " << aNorm << endl;
			
			//determine sign of a[1][1]
			bool aIsNegative = false;
			if (a.objMatrix[j][0] < 0) {
				aIsNegative = true;    
			}
			
			//compute the reflector
			if (aIsNegative) {
				v = a - e * aNorm;
			}
			else
				v = a + e * aNorm;
			
				//cout << "v= " << endl;
				//v.display();
			
			//compute v(Transpose) * v
			Matrix<T> vTransposev(1,1);
			Matrix<T> vTranspose(1,v.get_rows());
			vTranspose = v.transpose();
				//cout << "vtranspose= " << endl;
				//vtranspose.display();
			
			vTransposev = vTranspose * v;
			float vtv = vTransposev.objMatrix[0][0];
				//cout << "vtv= "<< vtv  << endl;
			
			Matrix temp(vTranspose.get_rows(), vTranspose.get_columns() );
			temp = vTranspose * (2/vtv);
			h = identity - ( v  * temp);
			
			
			Matrix<T> nextH(h.get_rows(), (*this).get_columns());
			nextH = h * (*this);
			(*this)= nextH;
			
			qRMatrix = h * qRMatrix;
			
			
		}
		
		
		(*this).replace_w_zeros();
		cout << "householder factorized matrix= " << endl;
		(*this).display();
		
		bColVector = qRMatrix * bColVector;
		cout << "bColVector after being multiplied by qRMatrix"<< endl;
		bColVector.display();
		
		Matrix householderFactorizedMatrix((*this).get_columns(),(*this).get_columns());
		for (int i = 0 ; i< (*this).get_columns(); i++) {
			for (int j =0 ; j < (*this).get_columns(); j++) {
				householderFactorizedMatrix.objMatrix[i][j] = (*this).objMatrix[i][j];
			}
		}
		
		Matrix<T> bColVecFinal(householderFactorizedMatrix.get_rows(),1);
		for (int i =0; i< householderFactorizedMatrix.get_rows(); i++) {
			bColVecFinal.objMatrix[i][0] = bColVector.objMatrix[i][0];
		}
		
		Matrix<T> answer(householderFactorizedMatrix.get_rows(),1);
		answer = householderFactorizedMatrix.back_substitution(bColVecFinal);
		
		cout << "shortened hh factorized matrix " << endl;
		householderFactorizedMatrix.display();
		
		cout << "bFinal = b shortened "<< endl;
		bColVecFinal.display();

		cout << "answer is "<< endl;
		answer.display();
		
		return answer;	
	}
	
	
	template <typename T>
	void Matrix<T>::householder()
	{
		
		int row = (*this).get_rows();
		int column = (*this).get_columns();
		(*this).display();
		Matrix<T> duplicate((*this).get_rows(), (*this).get_columns());
		duplicate = (*this);
		
		Matrix<T> qRMatrix(row,row);
		qRMatrix.set_identity();
		
		Matrix<T> h(row, row);
		for(int j=0; j<column  ; j++){
			//a e v h, same vectors/matricies from algorithm in notes
			Matrix<T> a(row, 1);
			Matrix<T> e(row, 1);
			Matrix<T> v(row, 1);
			
			
			//create identity matrix 
			Matrix<T> identity(row,row);
			identity.set_identity();
			
			
			int rowIndex = j;
			
			//make a = to current column
			for(int i=j; i<row; i++){
				a.objMatrix[i][0] = (*this).objMatrix[rowIndex++][j]; 
			}
			cout<< "a= " << endl;
			a.display();
			
			//initialize the e colomn vector
			e.objMatrix[j][0]=1;
			
			cout << "e= " << endl;
			e.display();
			
			//compute norm of current column
			float elementSqSum =0;
			for(int i=j; i<row; i++){
				elementSqSum += pow(a.objMatrix[i][0],2);
			}
			float aNorm = sqrt(elementSqSum);
			cout << "anorm= " << aNorm << endl;
			
			//determine sign of a11
			bool is_negative = false;
			if (a.objMatrix[j][0] < 0) {
				is_negative = true;    
			}
			
			//compute the reflector
			if (is_negative) {
				v = a - e * aNorm;
			}
			else
				v = a + e * aNorm;
			
			cout << "v= " << endl;
			v.display();
			
			//compute v(Transpose) * v
			Matrix<T> vtransposev(1,1);
			Matrix<T> vtranspose(1,v.get_rows());
			vtranspose = v.transpose();
			cout << "vtranspose= " << endl;
			vtranspose.display();
			
			vtransposev = vtranspose * v;
			float vtv = vtransposev.objMatrix[0][0];
			cout << "vtv= "<< vtv  << endl;
			
			Matrix temp(vtranspose.get_rows(), vtranspose.get_columns() );
			temp = vtranspose * (2/vtv);
			h = identity - ( v  * temp);
			
			cout <<"h= " << endl;
			h.replace_w_zeros();
			h.display();
			
			Matrix<T> h1a(h.get_rows(), (*this).get_columns());
			h1a = h * (*this);
			(*this)= h1a;
			
			cout << "h1a= " << endl;
			h1a.replace_w_zeros();
			h1a.display();
			
			
			cout << "this= "<< endl;
			(*this).display();
			
			qRMatrix = h * qRMatrix;
			
		}
		
		duplicate = qRMatrix * duplicate;
		
		cout << "duplicate is= "<< endl;
		duplicate.display();
		
		cout << "qRMatrix= " << endl;
		qRMatrix.display();
		
		cout << "H= " << endl;
		h.display();
		
		(*this).replace_w_zeros();
		cout << "ending this= " << endl;
		(*this).display();
		//return *this;	
		
		cout << "end" << endl;
		
	}//end householder
	
	
		
	template <typename T>
	T Matrix<T>::infinity_norm()	{
		
		T largestRowSum = 0;
		T currentRowSum = 0;
		for (int i =0; i< (*this).get_columns(); i++)
		{
			largestRowSum += abs((*this).objMatrix[0][i]);
		}
		
		for (int i =1; i< (*this).get_rows(); i++)
		{
			for (int j=0; j<(*this).get_columns(); j++) {
				currentRowSum += abs((*this).objMatrix[i][j]);
			}
			if (currentRowSum > largestRowSum) {
				largestRowSum = currentRowSum;
			}
			currentRowSum = 0;
		}
		
		return largestRowSum;
	
	}
	
	template <typename T>
	T Matrix<T>::get_max_i()	{
		T largestI = 0;
		T largestItem = abs((*this).objMatrix[0][0]);
		
		
		for (int i =1; i< (*this).get_rows(); i++)
		{
			if (abs((*this).objMatrix[i][0]) > largestItem) {
				largestItem = abs((*this).objMatrix[i][0]);
				largestI = i;
			}
		}
		
		return largestI;
		
	}
	
	template <typename T>
	T Matrix<T>::euclidean_norm()	
	{
		T sum =0;
		
		for (int i =0; i<(*this).rows_; i++) {
			sum += pow((*this).objMatrix[i][0], 2);
		}
		sum = sqrt(sum);
		return sum;
	
	}
	
	template <typename T>
	void Matrix<T>::replace_w_zeros()
	{ 
		for(int i=0; i< (*this).rows_; i++)
			for (int j=0; j< (*this).columns_; j++) 
			{
				if ((*this).objMatrix[i][j] < .001 && (*this).objMatrix[i][j] > -.001) 
				(*this).objMatrix[i][j] =0;
			}
	}
}//end namespace

#endif