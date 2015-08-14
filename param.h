/*
 *  param.h
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 5/7/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#include <vector>

#include "settings.h"
#include "path.h"

class MbRandom;
class path;
class settings;
class popsize;
class wfSamplePath;

class param {
	
public:
	param(MbRandom* r) {random = r; numProp = 0; numTunings = 0; numAccept = 0;};
	param(double x, MbRandom* r) {curVal = x; random = r; oldVal = x; tuning = 1; numTunings = 0; numProp = 0; numAccept = 0;}
	virtual double propose() = 0; //return proposal ratio!
	virtual double prior() = 0; //return prior ratio!
	virtual void updateTuning();
	void increaseProp() {numProp += 1;};
	void increaseAccept() {numAccept += 1;};
	virtual void reset() {curVal = oldVal;};
	double get() {return curVal;};
	double getOld() {return oldVal;};
	double getTuning() {return tuning;};
	void setOld(double v) {oldVal = v;};
	void setNew(double v) {curVal = v;};
	
	
protected:
	double curVal;
	double oldVal;
	double tuning;
	int numProp;
	int numAccept;
	int numTunings;
	MbRandom* random;
	
};

class param_gamma: public param {
public:
	param_gamma(double x, MbRandom* r): param(x, r) {scaling = 100.0; tuning=10.0;};
	double propose(); 
	double prior(); 
	
private:
	double scaling;
};

class param_h: public param {
public:
	param_h(double x, MbRandom* r): param(x, r) {scaling = 0.5; tuning=0.5;};
	double propose(); 
	double prior(); 
	
private:
	double scaling;
};

class param_path: public param {
public:
	param_path(path* p, param_gamma* al1, param_gamma* al2, MbRandom* r): param(r) {curPath = p; minUpdate = 10; fracOfPath = 10; min_dt = .001; grid = 10; a1 = al1; a2 = al2;};
	param_path(path* p, param_gamma* al1, param_gamma* al2, MbRandom* r, settings& s): param(r) {curPath = p; minUpdate = s.getMinUpdate(); fracOfPath = s.getFracOfPath(); min_dt = s.get_dt(); grid = s.get_grid(); fOrigin = acos(1.0-2.0*s.get_fOrigin()); a1 = al1; a2 = al2;};
	double propose();
	double proposeAlleleAge(double newAge);
	double proposeStart(double newStart);
	double proposeEnd(double newEnd);
	double propose(double x0, double xt, double t0, double t, std::vector<double> time_vec, int start_index, int end_index);
	double proposeAgePath(double x0,double xt,double t0,double t, std::vector<double> time_vec, int end_index);
	double prior() {return 0;};
	void updateTuning() {};
	void reset();
	path* get_path() {return curPath;};
	
private:
	int minUpdate;
	int fracOfPath;
	double min_dt;
	double grid;
	double fOrigin;
	path* curPath;
	path* newPath;
	path* oldPath;
	param_gamma* a1;
	param_gamma* a2;
	
	std::vector<double> make_time_vector(double newAge, int end_index, popsize* rho);
};

//NB: even though it says "freq", this really the transformed frequency
class start_freq: public param {
public:
	start_freq(double x, MbRandom* r, param_path* p): param(x, r) {curParamPath = p;};
	double propose();
	double prior(); 
	
private:
	param_path* curParamPath;
};

class end_freq: public param {
public:
	//SEE ABOVE!
	end_freq(double x, MbRandom* r, param_path* p): param(x, r) {curParamPath = p;};
	double propose();
	double prior();
	
private:
	param_path* curParamPath;
};

//if sample times are uncertain
class sample_time: public param {
	sample_time(double x, MbRandom* r): param(x, r) {};
	double propose();
	double prior();
};

//for learning the per-base-pair recombination rate
class param_rho: public param {
	param_rho(double x, MbRandom* r): param(x, r) {};
	double propose();
	double prior();
};

//for learning the age of the allele
class param_age: public param {
public:
	param_age(double x, MbRandom* r, param_path* p, double min_dt, double grid): param(x, r) {
		curParamPath = p; 
		popSize = ((wfSamplePath*)(p->get_path()))->get_pop();
		tuning = .2;
	};
	double propose();
	double prior();
	
private:
	param_path* curParamPath;
	popsize* popSize;
};