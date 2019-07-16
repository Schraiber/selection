/*
 *  measure.h
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 4/25/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#pragma once

#ifndef measure_H
#define measure_H

#include "math.h"
#include "gsl/gsl_sf_bessel.h"

#include "MbRandom.h"

class path;
class MbRandom;
class popsize;

class measure {
	
public:
	//constructor
	measure(MbRandom* r);
	
	//functions
	virtual double a(double x,double t) = 0; //drift function
	virtual double H(double x,double t) = 0; //potential function, int_0^x a(y,t)dy
	virtual double dadx(double x,double t) = 0; //derivative of drift function
	
	//path stuff
	//virtual path* prop_path(); //propose a path from the measure
	virtual path* prop_bridge(double x0, double xt, double t0, double t, std::vector<double>& time_vec) = 0; //propose a bridge
	//virtual void modify_path(path* p); //modify a path
	virtual double log_girsanov(path* p, measure* m, double lower, double upper, bool is_bridge = 0); //compute log likelihood ratio of path
	virtual double log_transition_density(double x, double y, double t) {return NAN;}; 
	
	//for bes0/wf stuff
	virtual double log_girsanov_wf(path* p, double alpha, bool is_bridge = 0) {return NAN;};
	virtual double log_girsanov_wf_r(path* p, double alpha, double h, popsize* rho, bool is_bridge = 0) {return NAN;};
	
protected:
	path* sim_bm(double x0, double t0, double t, std::vector<double>& time_vec);
	path* make_bb_from_bm(path* bm,double u, double v);
	MbRandom* random; 
	
};

//Measure of \arccos(1-2X_t) where X_t is a Wright-Fisher diffusion
class wfMeasure: public measure {
	
public:
	//constructor 
	wfMeasure(MbRandom* r, double g);
	//functions
	double a(double x, double t);
	double H(double x, double t);
	double dadx(double x, double t);
	//simulation
	path* prop_bridge(double x0, double xt, double t0, double t, std::vector<double>& time_vec, double rescale);
	path* prop_bridge(double x0, double xt, double t0, double t, std::vector<double>& time_vec) {return prop_bridge(x0, xt, t0, t,time_vec,-INFINITY);};
	//transform variable
	double fisher(double x) {return acos(1.0-2.0*x);};
	double inverse_fisher(double x) {return (1.0-cos(x))/2.0;};
	void invert_path(path* p); 
	//interact with selection coefficient
	double get_gamma() {return gamma;};
	void set_gamma(double g) {gamma = g;};
	//interact with things
	void set_num_test(int n) {num_test = n;};
	//allele age
	double expected_age(double f) {return -2.0*log(f)*f/(1.0-f); }; //returns the expected neutral age
	
private:
	double gamma; //selection coefficient
	int num_test; 
};

//Measure of 2\sqrt{X_t} where X_t is a critical continuous-state branching process
class cbpMeasure: public measure {
	
public:
	//constructor
	cbpMeasure(MbRandom* r): measure(r) {};
	//functions
	double a(double x, double t);
	double H(double x, double t);
	double dadx(double x, double t);
	//simulation
	path* prop_bridge(double x0, double xt, double t0, double t, std::vector<double>& time_vec);
	
	//transition density
	double log_transition_density(double x, double y, double t) {return log(x/t) - (x*x+y*y)/(2*t) + log(gsl_sf_bessel_I1_scaled(x*y/t))+x*y/t;};
	
	//for Wright-Fisher wrt Bessel measure with appropriate cancelling
	double a2_wf(double x, double t, double alpha); //NB: THIS IS ALREADY SQUARED!
	double H_wf(double x, double t, double alpha);
	double dadx_wf(double x, double t, double alpha);
	double log_girsanov_wf(path* p, double alpha, bool is_bridge = 0);
	
	//for Wright-Fisher relative to Wright-Fisher (in effect, Bessel measure cancelled out)
	//alpha1 is the new alpha, alpha2 is the old alpha
	double a2_wfwf(double x, double t, double alpha1, double alpha2);
	double H_wfwf(double x, double t, double alpha1, double alpha2);
	double dadx_wfwf(double x, double t, double alpha1, double alpha2);
	double log_girsanov_wfwf(path* p, double alpha1, double alpha2);
	
	//for Wright-Fisher variable population size wrt Bessel measure with appropriate cancelling
	double H_wf_r(double x, double t, double alpha, double h, popsize* rho, bool leftLimit = 0);
	double dHdt_wf_r(double x, double t, double alpha, double h, popsize* rho, bool leftLimit = 0);
	double a2_wf_r(double x, double t, double alpha, double h, popsize* rho, bool leftLimit = 0);
	double dadx_wf_r(double x, double t, double alpha, double h, popsize* rho, bool leftLimit = 0);
	double log_girsanov_wf_r(path* p, double alpha, double h, popsize* rho, bool is_bridge);
	
	//for Wright-Fisher with variable population size relative to Wright-Fisher with variable population size
	//alpha1 is new alpha, alpha 2 is old alpha
	double H_wfwf_r(double x, double t, double alpha1, double alpha1p, double alpha2, double alpha2p, popsize* rho, bool leftLimit = 0);
	double dHdt_wfwf_r(double x, double t, double alpha1, double alpha1p, double alpha2, double alpah2p, popsize* rho, bool leftLimit = 0);
	double a2_wfwf_r(double x, double t, double alpha1, double alpha1p, double alpha2, double alpha2p, popsize* rho, bool leftLimit = 0);
	double dadx_wfwf_r(double x, double t, double alpha1, double alpah1p, double alpha2, double alpha2p, popsize* rho, bool leftLimit = 0);
	double log_girsanov_wfwf_r(path* p, double alpha1, double alpha1p, double alpha2, double alpha2p, popsize* rho);

private:
	std::vector<double> rvMF(double kappa, int d); //generates a vonMises-Fisher random variable
	std::vector<double> unifSphere(int d); //generate a uniform random variable on the d-sphere
	double rW(double kappa, int m); //generate a random W, see Wood (1994)
};

//Measure of (\pi - 2\sqrt{X_t}) where X_t is a critical continuous-state branching process
class flippedCbpMeasure: public measure {

public:
	//constructor
	flippedCbpMeasure(MbRandom* r): measure(r) {};
	//functions
	double a(double x, double t);
	double H(double x, double t);
	double dadx(double x, double t);
	//simulation
	path* prop_bridge(double x0, double xt, double t0, double t, std::vector<double>& time_vec);
	
	//log transition density
	double log_transition_density(double x, double y, double t) {double w = PI - x; double z = PI - y; return log(w/t) - (w*w+z*z)/(2*t) + log(gsl_sf_bessel_I1_scaled(w*z/t))+w*z/t;};
	
	
};

class wienerMeasure: public measure {
	
public:
	//constructor
	wienerMeasure(MbRandom* r): measure(r) {};
	//functions 
	double a(double x, double t) {return 0;};
	double H(double x, double t) {return 0;};
	double dadx(double x, double t) {return 0;};
	
	//girsanov
	double log_girsanov(path* p, measure* m, double lower, double upper, bool is_bridge = 0);
	
	//transition density
	double log_transition_density(double x, double y, double t) {return -1.0/2.0*log(2*PI*t) - (x-y)*(x-y)/(2.0*t); };
	
	
	//simulate
	path* prop_bridge(double x0, double xt, double t0, double t, std::vector<double>& time_vec);
	path* prop_path(double x0, double t0, double t, std::vector<double>& time_vec);

private:
	path* make_bb_from_bm(path* bm,double u, double v);
	
};


#endif