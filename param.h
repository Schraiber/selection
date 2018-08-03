/*
 *  param.h
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 5/7/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#pragma once

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
    param(MbRandom* r) {curVal = 0; oldVal = 0; tuning = 1; random = r; numProp = 0; numTunings = 0; numAccept = 0; minTuning = 0;};
    param(double x, MbRandom* r) {curVal = x; random = r; oldVal = x; tuning = 1; numTunings = 0; numProp = 0; numAccept = 0; minTuning = 0;}
    ~param() {};
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
    void setNew(double v) {oldVal = curVal; curVal = v;};
	
	
protected:
	double curVal;
	double oldVal;
	double tuning;
    double minTuning;
	int numProp;
	int numAccept;
	int numTunings;
	MbRandom* random;
    
    //draw a reflected uniform RV
    //centered at x, of width w, but constrained to be between low and high
    double reflectedUniform(double x, double w, double low, double high);
	
};

class param_gamma: public param {
public:
    param_gamma(double x, MbRandom* r): param(x, r) {scaling = 10.0; tuning=10.0; minTuning=0.0;};
	double propose(); 
	double prior(); 
	
private:
	double scaling;
};

class param_h: public param {
public:
    param_h(double x, MbRandom* r): param(x, r) {scaling = 0.5; tuning=0.5; minTuning=0.0;};
	double propose(); 
	double prior(); 
	
private:
	double scaling;
};

class param_F: public param {
public:
    param_F(double x, MbRandom* r): param(x, r) {tuning = 0.5; minTuning = 0.0;};
    double propose();
    double prior();
    
};

class param_path: public param {
public:
	//param_path(path* p, param_gamma* al1, param_gamma* al2, MbRandom* r): param(r) {curPath = p; minUpdate = 10; fracOfPath = 10; min_dt = .001; grid = 10; a1 = al1; a2 = al2;};
	param_path(path* p, param_gamma* al1, param_gamma* al2, MbRandom* r, settings& s): param(r) {curPath = p; minUpdate = s.getMinUpdate(); fracOfPath = s.getFracOfPath(); min_dt = s.get_dt(); grid = s.get_grid(); fOrigin = acos(1.0-2.0*s.get_fOrigin()); a1 = al1; a2 = al2;};
    ~param_path() {delete curPath; delete newPath; delete oldPath;};
	double propose();
	double proposeAlleleAge(double newAge, double oldAge);
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
    start_freq(double x, MbRandom* r, param_path* p): param(x, r) {curParamPath = p; minTuning=0.0;};
	double propose();
	double prior();
    
    void reset();
	
private:
	param_path* curParamPath;
};

class end_freq: public param {
public:
	//SEE ABOVE!
    end_freq(double x, MbRandom* r, param_path* p): param(x, r) {curParamPath = p; minTuning=3.14/500.;};
	double propose();
	double prior();
    
    void reset();
	
private:
	param_path* curParamPath;
};

//if sample times are uncertain
class sample_time: public param {
public:
    sample_time(double x, double anc, double rec, int ss, int sc, MbRandom* r): param(x, r) {
        oldest = anc;
        youngest = rec;
        sample_size = ss;
        sample_count = sc;
        tuning = (youngest-oldest)/2.0;
        cur_idx = -1;
        old_idx = -1;
    };
    
    double get_ss() {return sample_size;};
    double get_sc() {return sample_count;};
    
    double get_oldest() {return oldest;};
    double get_youngest() {return youngest;};
    
    int get_oldest_idx() {return oldest_idx;};
    int get_youngest_idx() {return youngest_idx;};
    int get_idx() {return cur_idx;};
    
    void set_oldest_idx(int i) {oldest_idx = i;};
    void set_youngest_idx(int i) {youngest_idx = i;};
    void set_idx(int i) {old_idx = cur_idx; cur_idx = i;}; //keeps time the same, but changes the idx
    void reset_idx() {cur_idx = old_idx;};
    void reset();
    
    
	double propose();
	double prior();
    void updateTuning();
    
    void set_path(param_path* p) {curParamPath = p;};
    bool operator<(sample_time& t2) { return (curVal < t2.get()); };

private:
    //boundaries
    double oldest;
    double youngest;
    
    //index of most ancient time, most recent time, current value
    int oldest_idx;
    int youngest_idx;
    int cur_idx;
    int old_idx;
    
    //sample relevant information
    int sample_count;
    int sample_size;
    param_path* curParamPath;
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
        minTuning = min_dt*10;
	};
	double propose();
	double prior();
        
    void reset();
	
private:
	param_path* curParamPath;
	popsize* popSize;
};