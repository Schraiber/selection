/*
 *  settings.cpp
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 4/30/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#include "math.h"
#include "settings.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>

settings::settings(int argc, char* const argv[]) {
	
	//define defaults
	max_dt = .0001;
	min_grid = 10;
	bridge = 0;
	mcmc = 0;
	linked_sites = 0;
	num_gen = 100;
	p = 0;
	num_test = 1000;
	rescale = -INFINITY;
	printFreq = 1000;
	sampleFreq = 1000;
	fracOfPath = 20;
	min_update = 10;
	baseName = "out";
	fOrigin = 0.0;
	infer_age = 0;
	popFile = "";
	mySeed = time(0);
	output_tsv = 0;

	//read the parameters
	int ac = 1;
	while (ac < argc) {
		switch(argv[ac][1]) {
			case 'd':
				max_dt = atof(argv[ac+1]);
				ac += 2;
				break;
			case 'g':
				min_grid = atoi(argv[ac+1]);
				ac += 2;
				break;
			case 'n':
				num_gen = atoi(argv[ac+1]);
				ac += 2;
				break;
			case 'b':
				bridge = 1;
				bridge_pars = argv[ac+1];
				ac += 2;
				break;
			case 'R':
				output_tsv = 1;
				ac += 1;
				break;
			case 'p':
				p = 1;
				ac += 1;
				break;
			case 't':
				num_test = atoi(argv[ac+1]);
				ac += 2;
				break;
			case 'r':
				rescale = atof(argv[ac+1]);
				ac += 2;
				break;
			case 'X':
				mcmc = 1;
				sample_freq = argv[ac+1];
				ac += 2;
				break;
			case 'N':
				sample_size = argv[ac+1];
				ac += 2;
				break;
			case 'T':
				sample_time = argv[ac+1];
				ac += 2;
				break;
			case 'f':
				printFreq = atoi(argv[ac+1]);
				ac += 2;
				break;
			case 'o':
				baseName = std::string(argv[ac+1]);
				ac += 2;
				break;
			case 's':
				sampleFreq = atoi(argv[ac+1]);
				ac += 2;
				break;
			case 'F':
				fracOfPath = atoi(argv[ac+1]);
				ac += 2;
				break;
			case 'M':
				min_update = atoi(argv[ac+1]);
				ac += 2;
				break;
			case 'O':
				fOrigin = atof(argv[ac+1]);
				ac += 2;
				break;
			case 'e':
				mySeed = atoi(argv[ac+1]);
				ac += 2;
				break;
			case 'a':
				infer_age = 1;
				ac += 1;
				break;
			case 'P':
				popFile = argv[ac+1];
				ac += 2;
				break;
		}
	}
}

std::vector<double> settings::parse_bridge_pars() {
	std::vector<double> pars(0);
	std::string string_pars(bridge_pars);
	std::istringstream stringstream_pars(string_pars);
	std::string cur_par;
	
	while (std::getline(stringstream_pars,cur_par,',')) {
		pars.push_back(atof(cur_par.c_str()));
	}
	if (pars.size() < 4) {
		std::cout << "ERROR: Not enough bridge parameters" << std::endl;
		std::cout << "Only " << pars.size() << " specified; 4 are required" << std::endl;
		exit(1);
	}
	return pars;
}

void settings::print() {
	std::cout << "max_dt\t" << max_dt << std::endl;
	std::cout << "min_grid\t" << min_grid << std::endl;
	if (bridge) {
		std::cout << "num_test\t" << num_test << std::endl;
		std::cout << "rescale\t" << rescale << std::endl;
		std::vector<double> pars = parse_bridge_pars();
		std::cout << "x0\t" << pars[0] << std::endl;
		std::cout << "xt\t" << pars[1] << std::endl;
		std::cout << "gamma\t" << pars[2] << std::endl;
		std::cout << "t\t" << pars[3] << std::endl;
	} else if (mcmc) {
		std::cout << "num_gen\t" << num_gen << std::endl;
		if (linked_sites) {
			//file destinations
		} else {
			//parameters of the path
		}
	}
}
