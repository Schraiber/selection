/*
 *  matrix.cpp
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 4/26/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#include "matrix.h"

#include <iostream>

matrix::matrix(int nrow, int ncol) {
	double dat = 0;
	mat.resize(nrow);
	for (int i = 0; i < nrow; i++) {
		std::vector<double> cur_row(ncol,dat);
		mat[i] = cur_row;
	}
}

matrix::matrix(std::vector<double> vec, bool row) {
	if (row) {
		mat.resize(1);
		mat[0].resize(vec.size());
		for (int i = 0; i < vec.size(); i++) {
			mat[0][i] = vec[i];
		}
	} else {
		mat.resize(vec.size()); 
		for (int i = 0; i < vec.size(); i++) {
			mat[i][0] = vec[i];
		}
	}
}

matrix& matrix::operator=(matrix& b) {
	int nrow = b.get_nrow();
	mat.resize(nrow);
	for (int i = 0; i < nrow; i++) {
		mat[i] = b.row(i);
	}
	return *this;
}

matrix matrix::operator+(matrix& b) {
	int nrow = b.get_nrow();
	int ncol = b.get_ncol();
	if (mat.size() != nrow) {
		std::cout << "nrow(b) != my_nrow" << std::endl;
		exit(1);
	} else if (mat[0].size() != ncol) {
		std::cout << "ncol(b) != my_ncol" << std::endl;
		exit(1);
	}
	matrix new_matrix(nrow,ncol);
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			new_matrix(i,j) = mat[i][j] + b(i,j);
		}
	} 
	return new_matrix; 
}

matrix matrix::operator-(matrix& b) {
	int nrow = b.get_nrow();
	int ncol = b.get_ncol();
	if (mat.size() != nrow) {
		std::cout << "nrow(b) != mat.size()" << std::endl;
		exit(1);
	} else if (mat[0].size() != ncol) {
		std::cout << "ncol(b) != mat[0].size()" << std::endl;
		exit(1);
	}
	matrix new_matrix(nrow,ncol);
	for (int i = 0; i < nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			new_matrix(i,j) = mat[i][j]+b(i,j);
		}
	}
	return new_matrix;
}

matrix matrix::operator*(matrix& b) {
	int nrow = b.get_nrow();
	int ncol = b.get_ncol();
	if (mat[0].size() != nrow) {
		std::cout << "mat[0].size() != nrow, multiplication mismatch" << std::endl;
		exit(1);
	}
	matrix new_matrix(mat.size(),ncol);
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < ncol; j++) {
			for (int k = 0; k < nrow; k++) {
				new_matrix(i,j) += mat[i][k]*b(k,j);
			}
		}
	}
	return new_matrix;
}

matrix matrix::operator*(double& s) {
	matrix new_matrix(mat.size(), mat[0].size());
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			new_matrix(i,j) = mat[i][j]*s;
		}
	}
	return new_matrix;
}

matrix matrix::operator/(double& s) {
	matrix new_matrix(mat.size(), mat[0].size());
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			new_matrix(i,j) = mat[i][j]/s;
		}
	}
	return new_matrix;
}

matrix matrix::operator+(double& s) {
	matrix new_matrix(mat.size(), mat[0].size());
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			new_matrix(i,j) = mat[i][j]+s;
		}
	}
	return new_matrix;
}

matrix matrix::operator-(double& s) {
	matrix new_matrix(mat.size(), mat[0].size());
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			new_matrix(i,j) = mat[i][j]-s;
		}
	}
	return new_matrix;
}


matrix matrix::operator-() {
	matrix new_matrix(mat.size(), mat[0].size());
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			new_matrix(i,j) = -mat[i][j];
		}
	}
	return new_matrix;
}

void matrix::print(std::ostream& io) {
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			io << mat[i][j] << " ";
		}
		io << std::endl; 
	}
}

void matrix::t() {
	std::vector< std::vector< double > > tmp_mat = mat;
	for (int i = 0; i < mat.size(); i++) {
		for (int j = 0; j < mat[0].size(); j++) {
			mat[i][j] = tmp_mat[j][i];
		}
	}
}