/*
 *  matrix.h
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 4/26/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#ifndef matrix_H
#define matrix_H

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

class matrix {
	
public:
	//constructors
	matrix() {mat.resize(0);};
	matrix(int nrow, int ncol);
	matrix(std::vector<double> vec, bool row);
	
	//operators
	matrix& operator=(matrix& b); //copy
	matrix operator-(); //unary negative
	matrix operator+(matrix& b); //element-wise addition
	matrix operator+(double& s); //scalar addition
	matrix operator-(matrix& b); //element-wise subtraction
	matrix operator-(double& s); //scalar subtraction
	matrix operator*(matrix& b); //matrix multiplication
	matrix operator*(double& s); //scalar multiplication
	matrix operator/(double& s); //scalar divsion
	matrix operator%(matrix& b); //element-wise (hadamard) multiplication
	
	//elements
	double& operator()(int i, int j) {return mat[i][j];}; 
	std::vector<double>& row(int i) {return mat[i];};
	
	//operations
	void t();
	
	//information
	int get_nrow() {return mat.size();}; //gets number of rows
	int get_ncol() {return mat[0].size();}; //gets number of columns
	void print(std::ostream& io); 
	
	
private:
	std::vector< std::vector< double > > mat;
}; 



#endif