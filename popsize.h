/*
 *  popsize.h
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 11/22/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#ifndef popsize_H
#define popsize_H

#include <vector>
#include <string>
#include <stdlib.h>

class popsize {
public:
	//constructor takes in a path to a population size file
	popsize(std::string pop_size_file);
	
	//this gets the popsize at time t. Possibly the left limit
	double getSize(double t, bool leftLim = 0);
	//this gets the DERIVATIVE of the popsize at time t. Possibly the left limit.
	double getDeriv(double t, bool leftLim = 0);
	
	//this gets the transformed time. This is good
	double getTau(double t);
	std::vector<double> getTau(const std::vector<double>& t_vec);
	
	//gets the breakpoints spanned by an interval
	std::vector<double> getBreakTimes(double t0, double t);
	
	//get an entry of the vectors
	double getSizes(int i) {return sizes[i];};
	double getRates(int i) {return rates[i];};
	double getTimes(int i) {return times[i];};
	
	
private:
	//these functions all index the parameters of the piecewise exponential population history
	// nb: index 0 corresponds to the present, and all of the things are set to 0
	// nb: index n+1 corresponds to INFINITY, and rate should be 0
	// t[0] = 0 > t[1] > t[2] > ... > t[n] > t[n+1] = -INFINITY
	std::vector<double> sizes; //population size at most ancient time of an interval
	std::vector<double> rates; //growth rates in an interval
	std::vector<double> times; // begining times of intervals
	
	void computeT();
	std::vector<double> T; //these are the integrals over a whole interval
};

#endif