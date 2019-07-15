/*
 *  param.cpp
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 5/7/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#include "math.h"

#include "param.h"
#include "MbRandom.h"
#include "measure.h"
#include "path.h"
#include "popsize.h"

#include <algorithm>
#include <iomanip>

void param::updateTuning() {
	if (numProp > 0) {
		double curAccept = (double)numAccept/(double)numProp;
		double scaleFactor = 1.0/sqrt(numTunings);
		if (.01 < scaleFactor) {
			scaleFactor = .01;
		}
		if (curAccept > .3) {
			tuning = exp(log(tuning)+scaleFactor);
		} else {
			tuning = exp(log(tuning)-scaleFactor);
		}
        if (tuning < minTuning) {
            tuning = minTuning;
        }
	}
	numAccept = 0;
	numProp = 0;
	numTunings += 1;
}

double param::reflectedUniform(double x, double w, double low, double high) {
    double v;
    do {
        v = random->uniformRv(x-w/2.0,x+w/2.0);
        if (v < low) {
            v = 2*low-v;
        }
        else if (v > high) {
            v = 2*high-v;
        }
    } while (v < low || v > high);
    return v;
}

double param_gamma::propose() {
	oldVal = curVal;
	curVal = random->normalRv(oldVal,tuning);
	//double qOld = -1.0/2.0*log(2*PI*tuning*tuning) - (oldVal-curVal)*(oldVal-curVal)/(2.0*tuning*tuning);
	//double qNew = -1.0/2.0*log(2*PI*tuning*tuning) - (curVal-oldVal)*(curVal-oldVal)/(2.0*tuning*tuning);
	//return qOld-qNew;
    return 0;
}

double param_gamma::prior() {
	double pOld = -log(PI)+log(scaling)-log(oldVal*oldVal+scaling*scaling);
	double pNew = -log(PI)+log(scaling)-log(curVal*curVal+scaling*scaling);
	return pNew-pOld;
}

double param_h::propose() {
	oldVal = curVal;
	curVal = random->normalRv(oldVal,tuning);
	double qOld = -1.0/2.0*log(2*PI*tuning*tuning) - (oldVal-curVal)*(oldVal-curVal)/(2.0*tuning*tuning);
	double qNew = -1.0/2.0*log(2*PI*tuning*tuning) - (curVal-oldVal)*(curVal-oldVal)/(2.0*tuning*tuning);
	return qOld-qNew;
}

double param_h::prior() {
	double pOld = -log(PI)+log(scaling)-log((oldVal-0.5)*(oldVal-0.5)+scaling*scaling);
	double pNew = -log(PI)+log(scaling)-log((curVal-0.5)*(curVal-0.5)+scaling*scaling);
	return pNew-pOld;
}

double start_freq::propose() {
	oldVal = curVal;
    
    double propRatio = 0;
    
    //OLD: truncated normal
	//curVal = random->truncatedNormalRv(0, PI, oldVal, tuning);
	//double propRatio = random->truncatedNormalPdf(0, PI, curVal, tuning, oldVal);
	//propRatio -= random->truncatedNormalPdf(0, PI, oldVal, tuning, curVal);
    
    propRatio += curParamPath->proposeStart(curVal);

    //NEW: reflected uniform
    curVal = reflectedUniform(oldVal, tuning, 0, PI);
    propRatio += 0;
	
    return propRatio;
}

double start_freq::prior() {
	//uniform on [0,1] results in this density on the transformed space
	double pOld = log(sin(oldVal)) - log(2);
	double pNew = log(sin(curVal)) - log(2);
	return pNew - pOld;
	return 0;
}


double sample_time::propose() {
    //truncated normal
    oldVal = curVal;
    old_idx = cur_idx;
    
    //OLD: truncated normal
    curVal = random->truncatedNormalRv(oldest, youngest, oldVal, tuning);
    
    //NEW: reflected uniform
    //curVal = reflectedUniform(oldVal, tuning, oldest, youngest);
    
    double startVal = curVal;
    //Shift to closest value that's actually in the path
    //HOW BAD IS THIS IDEA???
    if (curVal < curParamPath->get_path()->get_time(0)) {
        //if it's older than the allele age, then set index = -1
        cur_idx = -1;
    } else if (curVal < oldVal) {
        //if less than, go down
        for (int i = old_idx; i >= 0; i--) {
            if (curParamPath->get_path()->get_time(i) < curVal) {
                double up_time = curParamPath->get_path()->get_time(i+1);
                double down_time = curParamPath->get_path()->get_time(i);
                double up_dif = up_time-curVal;
                double down_dif = curVal-down_time;
                //std::cout << "curVal = " << curVal << std::endl;
                //std::cout << "down_time = " << down_time << ", up_time = " << up_time << std::endl;
                //std::cout << "down_dif = " << down_dif << ", up_dif = " << up_dif << std::endl;
                if (down_dif < up_dif) {
                    curVal = down_time;
                    cur_idx = i;
                } else {
                    curVal = up_time;
                    cur_idx = i+1;
                }
                break;
            }
        }
    } else {
        //if greater than, go up
        for (int i = old_idx; i < curParamPath->get_path()->get_length(); i++) {
            if (i == -1) continue; //hack to deal with needing to start from the allele age
            if (curParamPath->get_path()->get_time(i) > curVal) {
                double up_time = curParamPath->get_path()->get_time(i);
                double down_time = curParamPath->get_path()->get_time(i-1);
                double up_dif = up_time-curVal;
                double down_dif = curVal-down_time;
                if (down_dif < up_dif) {
                    curVal = down_time;
                    cur_idx = i-1;
                } else {
                    curVal = up_time;
                    cur_idx = i;
                }
                break;
            }
        }
    }
    //try to be clever
//    if (cur_idx != -1 && sample_count > 0) {
//        ((wfSamplePath*)curParamPath->get_path())->updateFirstNonzero(curVal, oldVal);
//    }
    //just brute force...
    ((wfSamplePath*)curParamPath->get_path())->updateFirstNonzero();
    
    //OLD: truncated normal
    //double propRatio = random->truncatedNormalPdf(oldest, youngest, curVal, tuning, oldVal);
    //propRatio -= random->truncatedNormalPdf(oldest, youngest, oldVal, tuning, curVal);
    
    //NEW: refelcted uniform
    double propRatio = 0;
    if (curVal > youngest || curVal < oldest) {
        std::cout << "ERROR: sample_time proposal is outside of range" << std::endl;
        std::cout << "oldest = " << oldest << ", youngest = " << youngest << std::endl;
        std::cout << "Allele age = " << curParamPath->get_path()->get_time(0) << std::endl;
        std::cout << "oldVal = " << oldVal << ", curVal = " << curVal << std::endl;
        std::cout << "Starting curVal = " << startVal << std::endl;
        std::cout << "Final curVal = " << curVal << std::endl;
        std::cout << "old_idx = " << old_idx << ", cur_idx = " << cur_idx << std::endl;
        std::cout << "path->time(old_idx) = " << curParamPath->get_path()->get_time(old_idx) << std::endl;
        std::cout << "path->time(cur_idx) = " << curParamPath->get_path()->get_time(cur_idx) << std::endl;
        std::cout << "propRatio = " << propRatio << std::endl;
        //std::cin.ignore();
        //exit(1);
    }
    return propRatio;
}

double sample_time::prior() {
    return 0;
}

void sample_time::updateTuning() {
    param::updateTuning();
    if (tuning > youngest-oldest) {
        tuning = youngest-oldest;
    }
}

double param_age::propose() {
	oldVal = curVal;
	double topTime = ((wfSamplePath*)(curParamPath->get_path()))->get_firstNonzero();
    //OLD: truncated normal
	curVal = random->truncatedHalfNormalRv(topTime, 0, oldVal, tuning);
	double propRatio = log(random->truncatedHalfNormalPdf(topTime, 0, curVal, tuning, oldVal));
	propRatio -= log(random->truncatedHalfNormalPdf(topTime, 0, oldVal, tuning, curVal));
    //NEW: reflected uniform
    //curVal = reflectedUniform(oldVal, tuning, -INFINITY, topTime);
    //double propRatio = 0;
	if (propRatio != propRatio) {
		std::cout << "ERROR: Proposal ratio is nan! Debugging information:" << std::endl;
		std::cout << "oldVal: " << oldVal << " curVal: " << curVal << " tuning " << tuning << std::endl;
		std::cout << "log(P(theta | theta')) = " << log(random->truncatedHalfNormalPdf(topTime, 0, curVal, tuning, oldVal)) << std::endl;
		std::cout << "log(P(theta' | theta)) = " << log(random->truncatedHalfNormalPdf(topTime, 0, oldVal, tuning, curVal)) << std::endl;
	}
	propRatio += curParamPath->proposeAlleleAge(curVal, oldVal);
	return propRatio;
}

double param_age::prior() {
	//when popsize is bigger, have bigger likelihood of mutation
	double pOld = log(popSize->getSize(oldVal));
	double pNew = log(popSize->getSize(curVal));
	return pNew - pOld;
}

double end_freq::propose() {
	//truncated normal
	oldVal = curVal;
	curVal = random->truncatedNormalRv(0, PI, oldVal, tuning);
	double propRatio = random->truncatedNormalPdf(0, PI, curVal, tuning, oldVal);
	propRatio -= random->truncatedNormalPdf(0, PI, oldVal, tuning, curVal);
	propRatio += curParamPath->proposeEnd(curVal);
	return propRatio;
}

double end_freq::prior() {
	//no prior; implicit in the path!
	return 0;
}

//selects a random position to update
double param_path::propose() {	
	int start_index = random->discreteUniformRv(1, curPath->get_length()-(minUpdate+curPath->get_length()/fracOfPath));
	int end_index = start_index + minUpdate + curPath->get_length()/fracOfPath - 1; 
	double x0 = curPath->get_traj(start_index);
	double xt = curPath->get_traj(end_index);
	double t0 = curPath->get_time(start_index);
	double t = curPath->get_time(end_index);
	while (t - t0 < .0001 && end_index + minUpdate+curPath->get_length()/fracOfPath < curPath->get_length()) {
		end_index += minUpdate+curPath->get_length()/fracOfPath;
		t = curPath->get_time(end_index);
	}
	double propRatio = propose(x0,xt,t0,t,curPath->get_time(start_index,end_index),start_index,end_index);
	return propRatio;
}

//updates from the beginning
double param_path::proposeStart(double newStart) {
	int start_index = 0;
	int end_index = start_index + minUpdate+curPath->get_length()/fracOfPath;
	double x0 = newStart;
	double xt = curPath->get_traj(end_index);
	double t0 = curPath->get_time(start_index);
	double t = curPath->get_time(end_index);
	double propRatio = propose(x0,xt,t0,t,curPath->get_time(start_index,end_index),start_index,end_index);
	return propRatio;
}

//updates from the beginning
double param_path::proposeAlleleAge(double newAge, double oldAge) {
	int end_index;
	((wfSamplePath*)curPath)->set_update_begin();
    if (newAge < oldAge) {
        end_index = std::min(2*minUpdate, int(curPath->get_length()) - 1);
    } else {
        end_index = 0;
        while (curPath->get_time(end_index) < newAge) {
            end_index++;
        }
        if (end_index + 2*minUpdate < curPath->get_length()) {
            end_index += 2*minUpdate;
        }
        if (end_index > curPath->get_length()) {
            std::cout << "ERROR: trying to update allele age path past the end of the path!" << std::endl;
            std::cout << "path length = " << curPath->get_length() << ", end_index = " << end_index << std::endl;
            std::cout << "newAge = " << newAge << std::endl;
            curPath->print_time();
            exit(1);
        }
    }
	double x0 = fOrigin; 
	double t0 = newAge;
	double xt = curPath->get_traj(end_index);
	double t = curPath->get_time(end_index);
	while (t - t0 < 0.0001 && end_index + minUpdate+curPath->get_length()/fracOfPath < curPath->get_length()) {
		end_index += minUpdate+curPath->get_length()/fracOfPath;
		t = curPath->get_time(end_index);
	}
	popsize* rho = ((wfSamplePath*)curPath)->get_pop();
	std::vector<double> newTimeVector = make_time_vector(newAge, end_index, rho);
	double propRatio = proposeAgePath(x0,xt,t0,t,newTimeVector, end_index);
	return propRatio;
}

//this makes a time vector that hits the sample times and their boundaries
std::vector<double> param_path::make_time_vector(double newAge, int end_index, popsize* rho) {
	//figure out which times you need to include
	std::vector<double> timesToInclude;
	timesToInclude.push_back(newAge);
    
    //go through times to get which ones are in between
    double endTime = ((wfSamplePath*)curPath)->get_time(end_index);
    double curTime;
    double oldTime;
    double youngTime;
    sample_time* curSampleTime;
    for (int i = 0; i < ((wfSamplePath*)curPath)->get_num_samples(); i++) {
        curSampleTime = ((wfSamplePath*)curPath)->get_sampleTimeObj(i);
        oldTime = curSampleTime->get_oldest();
        curTime = curSampleTime->get();
        youngTime = curSampleTime->get_youngest();
        if (curTime >= newAge && curTime <= endTime) timesToInclude.push_back(curTime);
        if (oldTime >= newAge && oldTime <= endTime) timesToInclude.push_back(oldTime);
        if (youngTime >= newAge && youngTime <= endTime) timesToInclude.push_back(youngTime);
    }

    //generate the break times
    std::vector<double> breakPoints = rho->getBreakTimes(newAge,endTime);
    for (int j = 1; j < breakPoints.size(); j++) {
        timesToInclude.push_back(breakPoints[j]);
    }
    
    //sort
    std::sort(timesToInclude.begin(), timesToInclude.end());
    //ensure uniqueness
    std::vector<double>::iterator it = std::unique(timesToInclude.begin(), timesToInclude.end());
    timesToInclude.resize( std::distance(timesToInclude.begin(), it) );
    
    for (int j = 0; j < timesToInclude.size()-1; j++) {
        if (!(timesToInclude[j+1]>timesToInclude[j])) {
            std::cout << "ERROR: Times to include isn't strictly increasing!" << std::endl;
            for (int l = 0; l < timesToInclude.size(); l++) {
                std::cout << timesToInclude[l] << " ";
            }
            std::cout << std::endl;
            exit(1);
        }
    }
    
    
    //create the vector, going between each pair of things
	std::vector<double> newTimes;
    newTimes.push_back(timesToInclude[0]);
	for (int j = 0; j < timesToInclude.size()-1; j++) {
		double dt = min_dt;
		int steps = (timesToInclude[j+1]-timesToInclude[j])/dt+1;
        if (steps < minUpdate) {
            steps = minUpdate;
        }
        steps += 1;
        dt = (timesToInclude[j+1]-timesToInclude[j])/(steps-1);
        if (dt < std::numeric_limits<double>::epsilon()) {
            dt = std::numeric_limits<double>::epsilon();
            steps = (timesToInclude[j+1]-timesToInclude[j])/dt+1;
        }
        int end_k = newTimes.size()-1+steps;
        for (int k = newTimes.size(); k < end_k; k++) {
            newTimes.push_back(newTimes[k-1]+dt);
        }
        newTimes[newTimes.size()-1] = timesToInclude[j+1];
	}
    //check that time vector is strictly increasing
//    for (int j = 0; j < newTimes.size()-1; j++) {
//        if (!(newTimes[j+1]>newTimes[j])) {
//            std::cout << "ERROR: new time vector of length " << newTimes.size() << "  not strictly increasing" << std::endl;
//            std::cout << "Times to include are" << std::endl;
//            for (int l = 0; l < timesToInclude.size(); l++) {
//                std::cout << timesToInclude[l] << " ";
//            }
//            std::cout << std::endl;
//            std::cout << "newTimes[" << j << "+1] = " << newTimes[j+1] << ", newTimes[" << j << "] = " << newTimes[j] << std::endl;
//            exit(1);
//        }
//    }
	return newTimes;
}

//updates from the end
double param_path::proposeEnd(double newEnd) {
	int end_index = curPath->get_length()-1;
	int start_index = end_index - (minUpdate+curPath->get_length()/fracOfPath)+1;
	double x0 = curPath->get_traj(start_index);
	double xt = newEnd;
	double t0 = curPath->get_time(start_index);
	double t = curPath->get_time(end_index);
	double propRatio = propose(x0,xt,t0,t,curPath->get_time(start_index,end_index),start_index,end_index);
	return propRatio;
}

//does most of the hard work
double param_path::propose(double x0, double xt, double t0, double t, std::vector<double> time_vec, int start_index, int end_index) {
	//convert the times to tau times
	popsize* rho = ((wfSamplePath*)curPath)->get_pop();
	std::vector<double> tau_vec = rho->getTau(time_vec);
	double tau0 = rho->getTau(t0);
	double tau = rho->getTau(t);
	
	
	cbpMeasure myCBP(random);
	double dist_from_0 = x0;
	if (xt < x0) dist_from_0 = xt;
	double dist_from_pi = PI-xt;
	if (PI-x0 < PI-xt) dist_from_pi = PI-x0;
//	if (dist_from_0 < dist_from_pi) {
//		myCBP = new cbpMeasure(random);
//	} else {
//		myCBP = new flippedCbpMeasure(random);
//	}
	newPath = myCBP.prop_bridge(x0, xt, tau0, tau,tau_vec);
	oldPath = curPath->extract_path(start_index, end_index+1);
	
	newPath->replace_time(time_vec);
	curPath->modify(newPath,start_index);
	
	double propRatio = 0;
	
	//compute the likelihood ratio of current path under WF measure relative to CBP measure
	propRatio += myCBP.log_girsanov_wf_r(newPath, a1->get(), a2->get(),rho, 1);
	propRatio -= myCBP.log_girsanov_wf_r(oldPath, a1->get(), a2->get(),rho, 1);
	
	
	delete newPath;
	delete oldPath;
	
	return propRatio;
}

double param_path::proposeAgePath(double x0,double xt,double t0,double t, std::vector<double> time_vec, int end_index) {
	//convert the times to tau times
	popsize* rho = ((wfSamplePath*)curPath)->get_pop();
	std::vector<double> tau_vec = rho->getTau(time_vec);
	double tau0 = rho->getTau(t0);
	double tau = rho->getTau(t);
	
	cbpMeasure myCBP(random);
	newPath = myCBP.prop_bridge(x0, xt, tau0, tau, tau_vec);
	
	//these things, for computing the probability of the Bessel guy making it
	//should be in units of tau, so need to transform oldPath
	oldPath = curPath->extract_path(0,end_index+1);
	double tOld = rho->getTau(oldPath->get_time(oldPath->get_length()-1))-rho->getTau(oldPath->get_time(1));
	double tNew = newPath->get_time(newPath->get_length()-1)-newPath->get_time(1);
	
	newPath->replace_time(time_vec);
	((wfSamplePath*)curPath)->set_allele_age(t0, newPath, end_index);
	
	double propRatio = 0;
	
	//compute the likelihood ratio of current path under WF measure relative to CBP measure
	//NB: These ARE bridges but I want to compute the thing myself!
	propRatio += myCBP.log_girsanov_wf_r(newPath, a1->get(), a2->get(), rho,0);
	propRatio -= myCBP.log_girsanov_wf_r(oldPath, a1->get(), a2->get(), rho,0);
	
	propRatio += -1.0/2.0*xt*xt*(1.0/tNew-1.0/tOld)+2*log(tOld)-2*log(tNew);
	
	if (propRatio != propRatio) {
		std::cout << "ERROR: proposal ratio is nan! Debugging information sent to stderr:" << std::endl;
		std::cerr << "New path:" << std::endl;
		newPath->print_traj(std::cerr);
		newPath->print_time(std::cerr);
		std::cerr << myCBP.log_girsanov_wf_r(newPath, a1->get(), a2->get(), rho,0) << std::endl;
		std::cerr << "Old path:" << std::endl;
		oldPath->print_traj(std::cerr);
		oldPath->print_time(std::cerr);
		std::cerr << myCBP.log_girsanov_wf_r(oldPath, a1->get(), a2->get(), rho,0) << std::endl;
		std::cerr << "tNew tOld" << std::endl;
		std::cerr << tNew << " " << tOld << std::endl;
		std::cerr << "Time likelihood ratio" << std::endl;
		std::cerr << -1.0/2.0*xt*xt*(1.0/tNew-1.0/tOld)+2*log(tOld)-2*log(tNew) << std::endl;
		exit(1);
	}
	
	delete newPath;
	delete oldPath;
	
	return propRatio;
}

void param_path::reset() {
	((wfSamplePath*)curPath)->resetIntermediate();
}

void sample_time::reset() {
    curVal = oldVal;
    cur_idx = old_idx;
    ((wfSamplePath*)curParamPath->get_path())->resetFirstNonzero();
}

void param_age::reset() {
    curVal = oldVal;
    ((wfSamplePath*)curParamPath->get_path())->resetBeginning();
}

void end_freq::reset() {
    curVal = oldVal;
    ((wfSamplePath*)curParamPath->get_path())->resetIntermediate();
}

void start_freq::reset() {
    curVal = oldVal;
    ((wfSamplePath*)curParamPath->get_path())->resetIntermediate();
}