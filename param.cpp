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

double param_gamma::propose() {
	oldVal = curVal;
	curVal = random->normalRv(oldVal,tuning);
	double qOld = -1.0/2.0*log(2*PI*tuning*tuning) - (oldVal-curVal)*(oldVal-curVal)/(2.0*tuning*tuning);
	double qNew = -1.0/2.0*log(2*PI*tuning*tuning) - (curVal-oldVal)*(curVal-oldVal)/(2.0*tuning*tuning);
	return qOld-qNew;
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
	//truncated normal
	oldVal = curVal;
	curVal = random->truncatedNormalRv(0, PI, oldVal, tuning);
	double propRatio = random->truncatedNormalPdf(0, PI, curVal, tuning, oldVal);
	propRatio -= random->truncatedNormalPdf(0, PI, oldVal, tuning, curVal);
    propRatio += curParamPath->proposeStart(curVal);
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
    curVal = random->truncatedNormalRv(oldest, youngest, oldVal, tuning);
    double startVal = curVal;
    //Shift to closest value that's actually in the path
    //HOW BAD IS THIS IDEA???
    if (curVal < curParamPath->get_path()->get_time(0)) {
        //if it's older than the allele age, then set index = -1
        cur_idx = -1;
    } else if (curVal < oldVal) {
        //if less than, go down
        for (int i = old_idx; i > 0; i--) {
            if (curParamPath->get_path()->get_time(i) < curVal) {
                double up_time = curParamPath->get_path()->get_time(i+1);
                double down_time = curParamPath->get_path()->get_time(i);
                double up_dif = up_time-curVal;
                double down_dif = curVal-down_time;
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
    if (cur_idx!=-1) {
        ((wfSamplePath*)curParamPath->get_path())->updateFirstNonzero(curVal);
    }
    double propRatio = random->truncatedNormalPdf(oldest, youngest, curVal, tuning, oldVal);
    propRatio -= random->truncatedNormalPdf(oldest, youngest, oldVal, tuning, curVal);
    if (curVal > youngest || curVal < oldest) {
        std::cout << "Proposed new sample time" << std::endl;
        std::cout << "oldest = " << oldest << ", youngest = " << youngest << std::endl;
        std::cout << "Allele age = " << curParamPath->get_path()->get_time(0) << std::endl;
        std::cout << "oldVal = " << oldVal << ", curVal = " << curVal << std::endl;
        std::cout << "Starting curVal = " << startVal << std::endl;
        std::cout << "Final curVal = " << curVal << std::endl;
        std::cout << "propRatio = " << propRatio << std::endl;
        std::cin.ignore();
    }
    return propRatio;
}

double sample_time::prior() {
    return 0;
}

double param_age::propose() {
	oldVal = curVal;
	double topTime = ((wfSamplePath*)(curParamPath->get_path()))->get_firstNonzero();
	curVal = random->truncatedHalfNormalRv(topTime, 0, oldVal, tuning);
//  std::cout << "Vals: " << oldVal << " " << curVal << std::endl;
	double propRatio = log(random->truncatedHalfNormalPdf(topTime, 0, curVal, tuning, oldVal));
	propRatio -= log(random->truncatedHalfNormalPdf(topTime, 0, oldVal, tuning, curVal));
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
        end_index = 2*minUpdate;
    } else {
        end_index = 0;
        while (curPath->get_time(end_index) < newAge) {
            end_index++;
        }
        end_index += 2*minUpdate;
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

//this makes a time vector that hits the sample times
std::vector<double> param_path::make_time_vector(double newAge, int end_index, popsize* rho) {
	//figure out which times you need to include
	std::vector<double> timesToInclude;
	timesToInclude.push_back(newAge);
    
    //go through times to get which ones are in between
    double endTime = ((wfSamplePath*)curPath)->get_time(end_index);
    for (int i = 0; i < ((wfSamplePath*)curPath)->get_num_samples(); i++) {
        double curTime = ((wfSamplePath*)curPath)->get_sampleTimeValue(i);
        if (curTime >= newAge && curTime <= endTime) {
            timesToInclude.push_back(curTime);
        }
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
    
    //create the vector, going between each pair of things
	std::vector<double> newTimes;
	int cur_ind = 0;
	int temp_ind = 0;
	for (int j = 0; j < timesToInclude.size()-1; j++) {
		double dt = min_dt;
		int steps = (timesToInclude[j+1]-timesToInclude[j])/dt+1;
		newTimes.push_back(timesToInclude[j]);
		temp_ind++;
		for (int i = 1; i < steps - 1; i++) {
			newTimes.push_back(newTimes[cur_ind + i - 1] + dt);
			temp_ind++;
		}
		cur_ind += temp_ind;
		temp_ind = 0;
	}
	newTimes.push_back(timesToInclude[timesToInclude.size()-1]);
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
	
	
	measure* myCBP;
	double dist_from_0 = x0;
	if (xt < x0) dist_from_0 = xt;
	double dist_from_pi = PI-xt;
	if (PI-x0 < PI-xt) dist_from_pi = PI-x0;
	myCBP = new cbpMeasure(random);
//	if (dist_from_0 < dist_from_pi) {
//		myCBP = new cbpMeasure(random);
//	} else {
//		myCBP = new flippedCbpMeasure(random);
//	}
	newPath = myCBP->prop_bridge(x0, xt, tau0, tau,tau_vec);
	oldPath = curPath->extract_path(start_index, end_index+1);
	
	newPath->replace_time(time_vec);
	curPath->modify(newPath,start_index);
	
	double propRatio = 0;
	
	//compute the likelihood ratio of current path under WF measure relative to CBP measure
	propRatio += myCBP->log_girsanov_wf_r(newPath, a1->get(), a2->get(),rho, 1);
	propRatio -= myCBP->log_girsanov_wf_r(oldPath, a1->get(), a2->get(),rho, 1);
	
	
	delete myCBP;
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
	
	measure* myCBP;
	myCBP = new cbpMeasure(random); 
	newPath = myCBP->prop_bridge(x0, xt, tau0, tau, tau_vec);
	
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
	propRatio += myCBP->log_girsanov_wf_r(newPath, a1->get(), a2->get(), rho,0);
	propRatio -= myCBP->log_girsanov_wf_r(oldPath, a1->get(), a2->get(), rho,0);
	
	propRatio += -1.0/2.0*xt*xt*(1.0/tNew-1.0/tOld)+2*log(tOld)-2*log(tNew);
	
	if (propRatio != propRatio) {
		std::cout << "ERROR: proposal ratio is nan! Debugging information sent to stderr:" << std::endl;
		std::cerr << "New path:" << std::endl;
		newPath->print_traj(std::cerr);
		newPath->print_time(std::cerr);
		std::cerr << myCBP->log_girsanov_wf_r(newPath, a1->get(), a2->get(), rho,0) << std::endl;
		std::cerr << "Old path:" << std::endl;
		oldPath->print_traj(std::cerr);
		oldPath->print_time(std::cerr);
		std::cerr << myCBP->log_girsanov_wf_r(oldPath, a1->get(), a2->get(), rho,0) << std::endl;
		std::cerr << "Times New Old" << std::endl;
		std::cerr << tNew << " " << tOld << std::endl;
		std::cerr << "Time likelihood ratio" << std::endl;
		std::cerr << -1.0/2.0*xt*xt*(1.0/tNew-1.0/tOld)+2*log(tOld)-2*log(tNew) << std::endl;
		exit(1);
	}
	
	delete myCBP;
	delete newPath;
	delete oldPath;
	
	return propRatio;
}

void param_path::reset() {
	curPath->reset();
}

void sample_time::reset() {
    curVal = oldVal;
    cur_idx = old_idx;
    ((wfSamplePath*)curParamPath->get_path())->resetFirstNonzero();
}