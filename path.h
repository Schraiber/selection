/*
 *  path.h
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 4/25/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#pragma once

#ifndef path_H
#define path_H

#include <vector>
#include <string>
#include <iostream>
#include "math.h"

class settings;
class measure;
class wfMeasure;
class popsize;
class sample_time;
class MbRandom;

class path {

public:
	//constructor
	path() {trajectory.resize(0); time.resize(0); old_index = -1; old_trajectory.resize(0); old_time.resize(0);}; 
	path(double x0, double xt, double t0, double t, measure* m, settings& s);
	path(std::vector<double> p, std::vector<double> t) {trajectory = p; time = t;};
	path(double x0, double xt, double t0, double t, measure* m, std::vector<double> tvec);
	
	//element access
	std::vector<double> get_traj() {return trajectory;};
	double get_traj(int i) {return trajectory[i];};
	std::vector<double> get_traj(int i, int j);
	void set_traj(double x, int i) {trajectory[i] = x;};
	std::vector<double> get_time() {return time;};
	double get_time(int i) {return time[i];};
	std::vector<double> get_time(int i, int j);
	double get_length() {return trajectory.size();};
    double get_length_time() {return time.size();};
	std::vector<double>::iterator get_traj_iterator(int i) {return trajectory.begin()+i;};
	std::vector<double>::iterator get_time_iterator(int i) {return time.begin()+i;};
	path* extract_path(int i, int j);
	
	//element modification
	void flipCbp();
	void append(path* p); //adds the elements of p to the end of the current path
	void append(path* p, int i); //adds the elements of p starting with the ith element of p
	void insert(path* p, int i); //inserts the elements of p into the current path starting at index i of current path
	void modify(path* p, int i); //replaces current path with the elements of p starting at index i of current path
	virtual void reset(); //resets back to the stuff detailed in old_trajectory and old_time, starting from old_index
	void replace_time(std::vector<double> new_time); 
	
	//I/O
	void print(std::ostream& o = std::cout);
	void print_tsv(std::ostream& o = std::cout);
	virtual void print_traj(std::ostream& o = std::cout);
	void print_time(std::ostream& o = std::cout);
	
	void set_old_index(int i) {old_index = i;};
	
protected:
	//store the trajectory and the times
	//trajectory[i] corresponds to position at time[i]
	std::vector<double> trajectory;
	std::vector<double> time;	
	int old_index;
	std::vector<double> old_trajectory;
	std::vector<double> old_time;
};

//derived class that also has sample times and sample frequencies
//NB: sampleSizes and sampleCounts are in normal units!
class wfSamplePath : public path {
public:
	//constructor
    wfSamplePath(std::vector<double> p, std::vector<double> t) : path(p,t) {sample_time_vec.resize(0);};
	wfSamplePath(settings& s, wfMeasure* wf); //initializes a path from sample info, NB: does not propose the beginning!
    wfSamplePath(std::vector<sample_time*>& times, popsize* myPop, wfMeasure* wf, settings& s, MbRandom* r); //same as previous, but breaks out the parsing
	
	//destructor
	~wfSamplePath();
	
	//access sample aspects
	int get_num_samples() {return sample_time_vec.size();};
	int get_sampleTime(int i); //NOTE: RETURNS AN INDEX!
	double get_sampleSize(int i);
	double get_sampleCount(int i);
	double get_sampleFreq(int i);
	double get_firstNonzero();
    double get_sampleTimeValue(int i);
    sample_time* get_sampleTimeObj(int i);
    
    //change sample aspects
    void updateFirstNonzero(double t, double old_t);
    void updateFirstNonzero();
    void resetFirstNonzero() {first_nonzero=old_first_nonzero;};
	
	//for allele age stuff
	void set_allele_age(double a, path* p, int i); //this should set the allele age, prepend the new path starting at CURRENT i, and fix up sampleTime. 
	void set_update_begin(bool up = 1) {update_begin = up;}; //use this in the propose thing
	double get_allele_age() {return allele_age;};
	
	//reset the interior portion of the path
    void resetIntermediate();
    //reset the beginning of the path
    void resetBeginning();
    
	//sample probabilities
	double sampleProb(int i);
	std::vector<double> sampleProb();
    
    //ascertainment
    double ascertainModern(int min);
    double ascertainAncient();
    
	//print trajectory
	void print_traj(std::ostream& o = std::cout);
	
	//popsize
	popsize* get_pop() {return myPop;};
	
private:
	//relevant to the sample
    std::vector<sample_time*> sample_time_vec;
    std::vector<double> all_times;
	
	//relevant to allele age
	double allele_age; //the age itself
	double old_age; //the old age from the previous cycle
	bool update_begin; //was the begining modified?
	std::vector<double> old_begin_traj; //the trajectory from the old beginning
	std::vector<double> old_begin_time; //The times from the old beginning
	double first_nonzero; //the first sample time where there are more than 0 copies of the derived allele
    double old_first_nonzero;
	
	//parses comma separated list of parameters
	std::vector<double> parse_comma_sep(char* c);
    
    //parses input file
    void parse_input_file(std::string fin, int g, double N0);
    
    //for sorting input
    std::vector<int> orderTimeIndex();
    std::vector<double> sortByIndex(std::vector<double>& vec, std::vector<int> index);
    
	
	//the population size history
	popsize* myPop;

};

#endif 
