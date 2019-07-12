/*
 *  path.cpp
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 4/29/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#include "path.h"
#include "math.h"
#include "MbRandom.h"
#include "settings.h"
#include "measure.h"
#include "popsize.h"
#include "param.h"

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>

struct compare_index
{
    const std::vector<double> base_arr;
    compare_index (const std::vector<double> &arr) : base_arr (arr) {}
    
    bool operator () (int a, int b) const
    {
        return (base_arr[a] < base_arr[b]);
    }
};

//builds a bridge from x0 to xt 
path::path(double x0, double xt, double t0, double t, measure* m, settings& s) {
	//create the time vector
	double dt = s.get_dt();
	int min_steps = s.get_grid();
	int steps = (t-t0)/dt+1;
	if (steps < min_steps) {
		steps = min_steps;
	}
	steps += 1;
	dt = (t-t0)/(steps-1);
	time.resize(steps);
	time[0] = t0;
	for (int i = 1; i < steps; i++) {
		time[i] = time[i-1] + dt;
	}
	time[steps-1] = t; //HACK TO MAKE SURE THAT MACHINE ERROR DOESN'T FUCK ME UP
	path* temp = m->prop_bridge(x0, xt, t0, t,time);
	trajectory = temp->get_traj();
	delete temp;
}

//builds a bridge from x0 to xt with a fixed time vector
path::path(double x0, double xt, double t0, double t, measure* m, std::vector<double> tvec) {
	time = tvec;
	path* temp = m->prop_bridge(x0, xt, t0, t,time);
	trajectory = temp->get_traj();
	delete temp;
}

void path::print(std::ostream& o) {
	int i;
	o << "trajectory ";
	for (i = 0; i < trajectory.size(); i++) {
		o << trajectory[i] << " ";
	}
	o << std::endl;
	o << "time ";
	for (i = 0; i < trajectory.size(); i++) {
		o << time[i] << " ";
	}
	o << std::endl;
}

void path::print_tsv(std::ostream& o) {
	int i;
	o << "trajectory\ttime" << std::endl ;
	for (i = 0; i < trajectory.size(); i++) {
		o << trajectory[i] << "\t" << time[i] << std::endl;
	}
}

void path::print_traj(std::ostream& o) {
	int i;
	for (i = 0; i < trajectory.size(); i++) {
		o << trajectory[i] << " ";
	}
	o << std::endl;
}

void path::print_time(std::ostream& o) {
	int i;
	for (i = 0; i < time.size(); i++) {
		o << time[i] << " ";
	}
	o << std::endl;
}

void path::flipCbp() {
	for (int i = 0; i < trajectory.size(); i++) {
		trajectory[i] = PI-trajectory[i];
	}
}

void path::append(path* p) {
	for (int i = 0; i < p->get_length(); i++) {
		trajectory.push_back(p->get_traj(i));
		time.push_back(p->get_time(i));
	}
}

void path::append(path* p, int i) {
	for (int j = i; j < p->get_length(); j++) {
		trajectory.push_back(p->get_traj(j));
		time.push_back(p->get_time(j));
	}
}

void path::insert(path* p, int i) {
	trajectory.insert(trajectory.begin()+i,p->get_traj_iterator(0),p->get_traj_iterator(p->get_length()));
	time.insert(time.begin()+i,p->get_time_iterator(0),p->get_time_iterator(p->get_length()));
}

void path::modify(path* p, int i) {
	if (i != -1 || p != NULL) {
		old_index = i;
		int old_length = p->get_length();
		old_trajectory.resize(0);
		old_time.resize(0);
        //old_trajectory.resize(old_length);
        //old_time.resize(old_length);
		for (int j = 0; j < old_length; j++) {
			old_trajectory.push_back(trajectory[i+j]);
            //old_trajectory[j] = trajectory[i+j];
			trajectory[i+j] = p->get_traj(j);
			old_time.push_back(time[i+j]);
            //old_time[j] = time[i+j];
			time[i+j] = p->get_time(j);
            if (time[i+j] < time[i+j-1]) {
                std::cout << "ERROR: time vector is not sorted!" << std::endl;
                std::cout << time[i+j] << " >= " << time[i+j-1] << std::endl;
                exit(1);
            }
		}
        if (trajectory.size() != time.size()) {
            std::cout << "ERROR: Path trajectory and time are not same lenght!" << std::endl;
            std::cout << "trajectory.size() = " << trajectory.size() << std::endl;
            std::cout << "time.size() = " << time.size() << std::endl;
        }
	} else {
		old_trajectory.resize(0);
		old_time.resize(0);
		old_index = -1;
	}
}

void path::reset() {
	if (old_index != -1) {
		for (int j = 0; j < old_trajectory.size(); j++) {
			trajectory[old_index+j] = old_trajectory[j];
			time[old_index+j] = old_time[j];
		}
	}
    if (trajectory.size() != time.size()) {
        std::cout << "ERROR: Path trajectory and time are not same lenght!" << std::endl;
        std::cout << "trajectory.size() = " << trajectory.size() << std::endl;
        std::cout << "time.size() = " << time.size() << std::endl;
    }
}

std::vector<double> path::get_time(int i, int j) {
	std::vector<double> timeSlice(0,0);
	for (int k = i; k < j+1; k++) {
		timeSlice.push_back(time[k]);
	}
	return timeSlice;
}

std::vector<double> path::get_traj(int i, int j) {
	std::vector<double> trajSlice(0,0);
	for (int k = i; k < j+1; k++) {
		trajSlice.push_back(trajectory[k]);
	}
	return trajSlice;
}

path* path::extract_path(int i, int j) {
	std::vector<double> new_time(0);
	std::vector<double> new_traj(0);
	for (int k = i; k < j; k++) {
		new_time.push_back(time[k]);
		new_traj.push_back(trajectory[k]);
	}
	path* new_path = new path(new_traj,new_time);
	return new_path;
}

std::vector<double> wfSamplePath::parse_comma_sep(char* c) {
	std::vector<double> pars(0);
	std::string string_pars(c);
	std::istringstream stringstream_pars(string_pars);
	std::string cur_par;
	
	while (std::getline(stringstream_pars,cur_par,',')) {
		pars.push_back(atof(cur_par.c_str()));
	}
	return pars;
}


std::vector<double> wfSamplePath::sortByIndex(std::vector<double>& vec, std::vector<int> index) {
    std::vector<double> temp_vec(index.size());
    for (int i = 0; i < index.size(); i++) {
        temp_vec[i] = vec[index[i]];
    }
    return temp_vec;
}

double wfSamplePath::sampleProb(int i) {
	int idx = sample_time_vec[i]->get_idx();
    double sc = sample_time_vec[i]->get_sc();
    double ss = sample_time_vec[i]->get_ss();
	double sp = 0;
	if (idx != -1 && idx != 0) {
		sp += lgamma(ss+1)-lgamma(sc+1)-lgamma(ss-sc+1);
		sp += sc*log((1.0-cos(trajectory[idx]))/2.0);
		sp += (ss-sc)*log(1-(1.0-cos(trajectory[idx]))/2.0);
	} else {
		if (sc == 0) {
			sp += 0;
		} else {
			sp += -INFINITY;
		}
	}
	return sp;
}

std::vector<double> wfSamplePath::sampleProb() {
	std::vector<double> sp(0);
	for (int i = 0; i < sample_time_vec.size(); i++) {
		sp.push_back(sampleProb(i));
	}
	return sp;
}

void wfSamplePath::print_traj(std::ostream& o) {
	int i;
	for (i = 0; i < trajectory.size(); i++) {
		o << (1.0-cos(trajectory[i]))/2.0 << " ";
	}
	o << std::endl;
}

wfSamplePath::wfSamplePath(std::vector<sample_time*>& st, popsize* p, wfMeasure* wf, settings& s, MbRandom* r): path() {

    std::cout << "Creating initial path" << std::endl;
    
    myPop = p;
    
    sample_time_vec = st;
    
    int num_samples = st.size();
    
    //Initialize path
    std::vector<double> initial_data(num_samples);
    first_nonzero = -INFINITY;
    for (int i = 0; i < num_samples; i++) {
        //draw initial state from beta distribution with prior biased toward low frequency
        initial_data[i] = r->betaRv(.5 + sample_time_vec[i]->get_sc(), 2 + sample_time_vec[i]->get_ss() - sample_time_vec[i]->get_sc());
        if (first_nonzero == -INFINITY && sample_time_vec[i]->get_sc() != 0) {
            first_nonzero = sample_time_vec[i]->get();
        }
    }
    
    std::cout << "First nonzero timepoint is " << first_nonzero << std::endl;
    
    //get the times to include
    std::vector<double> breakPoints;
    for (int i = 0; i < num_samples; i++) {
        breakPoints.push_back(sample_time_vec[i]->get_youngest());
        breakPoints.push_back(sample_time_vec[i]->get());
        breakPoints.push_back(sample_time_vec[i]->get_oldest());
    }

    double min_time = *std::min_element(breakPoints.begin(), breakPoints.end());
    double max_time = *std::max_element(breakPoints.begin(), breakPoints.end());

    std::vector<double> curBreaks = myPop->getBreakTimes(std::min(min_time, first_nonzero), max_time);

    for (int i = 0; i < curBreaks.size(); i++) {
        breakPoints.push_back(curBreaks[i]);
    }
    
    //sort
    std::sort(breakPoints.begin(), breakPoints.end());
    //ensure uniqueness
    std::vector<double>::iterator it = std::unique(breakPoints.begin(), breakPoints.end());
    breakPoints.resize( std::distance(breakPoints.begin(), it) );
    
    //create the time vector
    double dt = s.get_dt();
    int min_steps = s.get_grid();
    int cur_end_ind = 0;
    int curBreakStart = 0;
    path* nextPath;
    
    
    int cur_time_idx = 1;
    int startBreak = 0;
    if (breakPoints[0] != sample_time_vec[0]->get()) {
        startBreak = 1;
    }
    std::vector<double> time_vec;
    time_vec.resize(0);
    time_vec.push_back(breakPoints[startBreak]);
    double curStart;
    double curEnd;
    for (int curBreak = startBreak; curBreak < breakPoints.size()-1; curBreak++) {
        //make sure the time vector includes all the break points
        curStart = breakPoints[curBreak];
        curEnd = breakPoints[curBreak+1];
        int steps = (curEnd-curStart)/dt+1;
        if (steps < min_steps) {
            steps = min_steps;
        }
        steps += 1;
        dt = (curEnd-curStart)/(steps-1);
        if (dt < std::numeric_limits<double>::epsilon()) {
            dt = std::numeric_limits<double>::epsilon();
            steps = (curEnd-curStart)/dt+1;
        }
        cur_end_ind++;
        int end_k = time_vec.size()-1+steps;
        for (int k = time_vec.size(); k < end_k; k++) {
            time_vec.push_back(time_vec[k-1]+dt);
            cur_end_ind++;
        }
        time_vec[time_vec.size()-1] = curEnd;
        //if we hit the time of a data point, simulate path between the two data points
        if (curEnd == sample_time_vec[cur_time_idx]->get()) {
            nextPath = new path(wf->fisher(initial_data[curBreakStart]), wf->fisher(initial_data[curBreakStart+1]), time_vec[0], time_vec[time_vec.size()-1], wf, time_vec);
            if (curBreakStart == 0) {
                this->append(nextPath);
            } else {
                this->append(nextPath,1);
            }
            sample_time_vec[cur_time_idx]->set_idx(cur_end_ind-1);
            cur_time_idx++;
            curBreakStart++;
            delete nextPath;
            time_vec.resize(0);
            time_vec.push_back(curEnd);
            
        } 
        cur_end_ind--;
    }
    
    std::cout << "Finished creating initial path" << std::endl;
}


wfSamplePath::~wfSamplePath() {
	delete myPop;
}

void wfSamplePath::set_allele_age(double a, path* p, int i) {
    int oldLength = time.size();
    double endTimeUpdate = p->get_time(p->get_length()-1);
	old_age = allele_age;
	allele_age = a;
	old_begin_traj = this->get_traj(0,i);
	old_begin_time = this->get_time(0,i);
	std::vector<double> tempTraj = p->get_traj();
	std::vector<double> tempTime = p->get_time();
	old_index = tempTime.size() - 1; 
	for (int j = i+1; j < trajectory.size(); j++) {
		tempTraj.push_back(trajectory[j]);
		tempTime.push_back(time[j]);
	}
	trajectory = tempTraj;
	time = tempTime;
    int newLength = time.size();
    int lengthDif = newLength - oldLength;
    
    //go through sample times to update index
    std::vector<double>::iterator search_it;
    int new_idx;
    for (int j = 0; j < sample_time_vec.size(); j++) {
        if (sample_time_vec[j]->get() < allele_age) {
            //if it's older than the allele age, just set idx to -1
            sample_time_vec[j]->set_idx(-1);
        } else if (sample_time_vec[j]->get() <= endTimeUpdate) {
            //if it's part of the new trajectory, find where it is
            search_it = std::lower_bound(time.begin(),time.end(),sample_time_vec[j]->get());
            if (search_it == time.end() && sample_time_vec[j]->get() != time[time.size()-1]) {
                std::cout << "ERROR: could not find sample time index " << j << " with value " << sample_time_vec[j]->get() << " in time vector!" << std::endl;
            }
            new_idx = search_it-time.begin();
            sample_time_vec[j]->set_idx(new_idx);
        } else {
            //otherwise, just update the position by adding!
            new_idx = sample_time_vec[j]->get_idx() + lengthDif;
            sample_time_vec[j]->set_idx(new_idx);
        }

    }
}



void wfSamplePath::resetIntermediate() {
    //check some things
    if (update_begin) {
        std::cout << std::endl << "ERROR: Trying to reset an intermediate part of the path, but allele age was updated!" << std::endl;
        exit(1);
    }
    if (old_index == -1) {
        std::cout << std::endl << "ERROR: Trying to reset an intermdiate part of the path, but old_index = -1!" << std::endl;
        exit(1);
    }
    //replace the trajectory
    for (int j = 0; j < old_trajectory.size(); j++) {
        trajectory[old_index+j] = old_trajectory[j];
        time[old_index+j] = old_time[j];
    }
    old_index = -1;
}

void wfSamplePath::resetBeginning() {
    //check some things
    if (!update_begin) {
        std::cout << std::endl << "ERROR: Trying to reset the beginning of the path, but allele age wasn't updated!" << std::endl;
        exit(1);
    }
    
    //prepend the old begining to the rest of the path
    allele_age = old_age;
    std::vector<double> tempTraj = old_begin_traj;
    std::vector<double> tempTime = old_begin_time;
    for (int j = old_index+1; j < trajectory.size(); j++) {
        tempTraj.push_back(trajectory[j]);
        tempTime.push_back(time[j]);
    }
    trajectory = tempTraj;
    time = tempTime;
        
    //also reset all the indices of the sample times
    for (int i = 0; i < sample_time_vec.size(); i++) {
        sample_time_vec[i]->reset_idx();
    }
    
    update_begin = 0;
}

int wfSamplePath::get_sampleTime(int i) {
    return sample_time_vec[i]->get_idx();
}

double wfSamplePath::get_sampleSize(int i) {
    return sample_time_vec[i]->get_ss();
}

double wfSamplePath::get_sampleCount(int i) {
    return sample_time_vec[i]->get_sc();
}

double wfSamplePath::get_sampleFreq(int i) {
    return (1.0-cos(trajectory[sample_time_vec[i]->get_idx()]))/2.0;
}

double wfSamplePath::get_firstNonzero() {
    return first_nonzero;
}

double wfSamplePath::get_sampleTimeValue(int i) {
    return sample_time_vec[i]->get();
}

sample_time* wfSamplePath::get_sampleTimeObj(int i) {
    return sample_time_vec[i];
}

void wfSamplePath::updateFirstNonzero(double t, double old_t) {
    old_first_nonzero = first_nonzero;
    if (t < old_first_nonzero || old_first_nonzero == old_t) {
        first_nonzero = t;
    }
}

void wfSamplePath::updateFirstNonzero() {
    old_first_nonzero = first_nonzero;
    first_nonzero = 0;
    for (int i = 0; i < sample_time_vec.size(); i++) {
        if (sample_time_vec[i]->get_sc() > 0 && sample_time_vec[i]->get() < first_nonzero) {
            first_nonzero = sample_time_vec[i]->get();
        } else {
        }
    }
}

void path::replace_time(std::vector<double> new_time) {
	if (new_time.size() != time.size()) {
		std::cout << "ERROR: Trying to replace a time vector with one of a different size!" << std::endl;
		exit(1);
	}
	time = new_time;
}
