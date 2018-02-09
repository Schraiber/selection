/*
 *  settings.h
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 4/30/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#pragma once

#ifndef settings_H
#define settings_H

#include <string>
#include <vector>
#include <stdlib.h>

class sample_time;
class MbRandom;
class popsize;

class settings {
public:
	settings(int argc, char* const argv[]);
	
	//print settings
	void print();
	
	//get stuff
	double get_dt() {return max_dt;};
	int get_grid() {return min_grid;};
	bool get_bridge() {return bridge;};
	bool get_mcmc() {return mcmc;};
	bool get_linked() {return linked_sites;};
	int get_num_gen() {return num_gen;};
	bool get_p() {return p;};
	int get_num_test() {return num_test;};
	double get_rescale() {return rescale;};
	int get_printFreq() {return printFreq;};
	int get_sampleFreq() {return sampleFreq;};
	std::string get_baseName() {return baseName;};
	int getMinUpdate() {return min_update;};
	int getFracOfPath() {return fracOfPath;};
	double get_fOrigin() {return fOrigin;};
	double get_seed() {return mySeed;};
	bool get_infer_age() {return infer_age;};
	std::string get_popFile() {return popFile;};
	bool get_output_tsv() {return output_tsv;};
    double get_a1prop() {return a1prop;};
    double get_a2prop() {return a2prop;};
    double get_ageprop() {return ageprop;};
    double get_endprop() {return endprop;};
    double get_timeprop() {return timeprop;};
    double get_pathprop() {return pathprop;};
    double get_a1start() {return a1start;};
    double get_a2start() {return a2start;};
    bool get_set_gen() {return set_gen;};
    bool get_set_N0() {return set_N0;};
    int get_gen_time() {return gen_time;};
    double get_N0() {return N0;};
    std::string get_infile() {return inputFile;};
		
	//parse things
	std::vector<double> parse_bridge_pars();
    std::vector<sample_time*> parse_input_file(MbRandom* r);
    popsize* parse_popsize_file();
	
private:
	double max_dt; //largest acceptable dt for paths
	int min_grid; //smallest acceptable number of grid points between the start and end of a path
	int num_gen; //number of mcmc cyles to run
	bool mcmc; //do mcmc?
	bool linked_sites; //incorproate linked sites?
	bool bridge; //rejection sample a bridge?
	bool p; //print the settings to stdout?
	char* bridge_pars; //parameters for diffusion bridge: x0,xt,gamma,t
	int num_test; //number of tests to try to estimate the scaling factor for rejection sampling
	double rescale; //scaling factor for rejection sampling
	int printFreq;
	std::string baseName;
	int sampleFreq; //how often to sample the chain. NOTE DIFFERENCE FROM sample_freq!!!
	int fracOfPath;
	int min_update;
	double fOrigin; //the frequency of a new mutation---should not matter!
	double mySeed; //for initializing mbRandom.
	bool infer_age; //should we infer the allele age
	std::string popFile;
	int output_tsv; //whether to output -b in tsv format
    double a1prop;
    double a2prop;
    double ageprop;
    double endprop;
    double pathprop;
    double timeprop;
    double a1start;
    double a2start;
    bool set_gen;
    bool set_N0;
    int gen_time;
    double N0;
    std::string inputFile;
};


#endif
