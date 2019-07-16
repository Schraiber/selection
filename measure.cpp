/*
 *  measure.cpp
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 4/25/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#include <vector>
#include <cmath>
#include<iomanip>
#include <signal.h>

#include "measure.h"
#include "path.h"
#include "MbRandom.h"
#include "popsize.h"

measure::measure(MbRandom* r) {
	random = r;
}

wfMeasure::wfMeasure(MbRandom* r, double g) : measure(r) {
	gamma = g;
	num_test = 0;
}

//uses Ito's formula to avoid computing an Ito integral
double measure::log_girsanov(path* p, measure* m,double lower, double upper, bool is_bridge) {
	int i;
	int path_len = p->get_length();
	double x0 = p->get_traj(0);
	double xt = p->get_traj(path_len - 1);
	double t0 = p->get_time(0);
	double tt = p->get_time(path_len -1); 
	
	//check if the path is in the bounds
	for (i = 0; i < path_len; i++) {
		if (p->get_traj(i) < lower || p->get_traj(i) > upper) {
			return -INFINITY;
		}
	}
	
	//compute everything for current measure (dominating measure)
	double H_wt = H(xt,tt);
	double H_w0 = H(x0,t0);
	//compute the time integral of the derivative
	double int_deriv = 0;
	for (i = 1; i < path_len; i++) {
		//use trapezoid rule: average the values at t_i, t_{i-1}
		int_deriv += (dadx(p->get_traj(i),p->get_time(i))+dadx(p->get_traj(i-1),p->get_time(i-1)))/2.0
					* (p->get_time(i)-p->get_time(i-1));
	}
	//compute the time integral of the square
	double int_square = 0;
	for (i = 1; i < path_len; i++) {
		int_square += (pow(a(p->get_traj(i),p->get_time(i)),2)+pow(a(p->get_traj(i-1),p->get_time(i-1)),2))/2.0
		              * (p->get_time(i)-p->get_time(i-1));
	}
	
	//compute everything for test measure m
	double Hm_wt = m->H(p->get_traj(path_len-1),p->get_time(path_len-1));
	double Hm_w0 = m->H(p->get_traj(0), p->get_time(0));
	//compute the time integral of the derivative
	double int_mderiv = 0;
	for (i = 1; i < path_len; i++) {
		//use trapezoid rule: average the values at t_i, t_{i-1}
		int_mderiv += (m->dadx(p->get_traj(i),p->get_time(i))+m->dadx(p->get_traj(i-1),p->get_time(i-1)))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}
	//compute the time integral of the square
	double int_msquare = 0;
	for (i = 1; i < path_len; i++) {
		int_msquare += (pow(m->a(p->get_traj(i),p->get_time(i)),2)+pow(m->a(p->get_traj(i-1),p->get_time(i-1)),2))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}
	
	if (!is_bridge) {
		return (Hm_wt-Hm_w0-1.0/2.0*int_mderiv-1.0/2.0*int_msquare)-(H_wt-H_w0-1.0/2.0*int_deriv-1.0/2.0*int_square);
	} else {
		double gir = (Hm_wt-Hm_w0-1.0/2.0*int_mderiv-1.0/2.0*int_msquare)-(H_wt-H_w0-1.0/2.0*int_deriv-1.0/2.0*int_square);
		//Computes the factors necessary for conditioning. NB: dominating measure ON TOP.
		double cond = log_transition_density(x0, xt, tt-t0)-m->log_transition_density(x0, xt, tt-t0);
		return gir + cond;
	}
}

//log girsanov for wiener measure; don't waste time computing the stuff for the wiener measure
double wienerMeasure::log_girsanov(path* p, measure* m, double lower, double upper, bool is_bridge) {
	int i;
	int path_len = p->get_length();
	double x0 = p->get_traj(0);
	double xt = p->get_traj(path_len - 1);
	double t0 = p->get_time(0);
	double tt = p->get_time(path_len -1);
		
	//compute everything for test measure m
	double Hm_wt = m->H(xt,tt);
	double Hm_w0 = m->H(x0, t0);
	//compute the time integral of the derivative
	double int_mderiv = 0;
	for (i = 1; i < path_len; i++) {
		//also check if legal
		if (p->get_traj(i) < lower || p->get_traj(i) > upper) {
			return -INFINITY;
		}
		//use trapezoid rule: average the values at t_i, t_{i-1}
		int_mderiv += (m->dadx(p->get_traj(i),p->get_time(i))+m->dadx(p->get_traj(i-1),p->get_time(i-1)))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}
	//compute the time integral of the square
	double int_msquare = 0;
	for (i = 1; i < path_len; i++) {
		int_msquare += (pow(m->a(p->get_traj(i),p->get_time(i)),2)+pow(m->a(p->get_traj(i-1),p->get_time(i-1)),2))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}
		
	if (!is_bridge) {
		return (Hm_wt-Hm_w0-1.0/2.0*int_mderiv-1.0/2.0*int_msquare);
	} else {
		double gir = (Hm_wt-Hm_w0-1.0/2.0*int_mderiv-1.0/2.0*int_msquare);
		double cond = -m->log_transition_density(x0, xt, tt-t0) + log_transition_density(x0, xt, tt-t0);
		return gir + cond;
	}
}

double cbpMeasure::log_girsanov_wf(path* p, double alpha, bool is_bridge) {
	int i;
	int path_len = p->get_length();
	double x0 = p->get_traj(0);
	double xt = p->get_traj(path_len - 1);
	double t0 = p->get_time(0);
	double tt = p->get_time(path_len -1);
	
	//compute everything for test measure m
	double Hm_wt = H_wf(xt,tt,alpha);
	double Hm_w0 = H_wf(x0, t0,alpha);
	//compute the time integral of the derivative
	double int_mderiv = 0;
	for (i = 1; i < path_len; i++) {
		//use trapezoid rule: average the values at t_i, t_{i-1}
		int_mderiv += (dadx_wf(p->get_traj(i),p->get_time(i),alpha)+dadx_wf(p->get_traj(i-1),p->get_time(i-1),alpha))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}
	//compute the time integral of the square
	double int_msquare = 0;
	for (i = 1; i < path_len; i++) {
		int_msquare += (a2_wf(p->get_traj(i),p->get_time(i),alpha)+a2_wf(p->get_traj(i-1),p->get_time(i-1),alpha))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}
	
	if (!is_bridge) {
		return (Hm_wt-Hm_w0-1.0/2.0*int_mderiv-1.0/2.0*int_msquare);
	} else {
		double gir = (Hm_wt-Hm_w0-1.0/2.0*int_mderiv-1.0/2.0*int_msquare);
		double cond = log_transition_density(x0, xt, tt-t0);
		return gir + cond;
	}
}

double cbpMeasure::log_girsanov_wf_r(path* p, double alpha1, double alpha2, popsize* rho, bool is_bridge) {
	int i;
	int j;
	int path_len = p->get_length();
	double x0 = p->get_traj(0);
	double xt = p->get_traj(path_len - 1);
	//these are in regular time
	double t0 = p->get_time(0);
	double tt = p->get_time(path_len -1);
	//convert to tau
	double tau0 = rho->getTau(t0);
	double taut = rho->getTau(tt);
	
	//compute everything for test measure m
	double Hm_wt = 0;
	double Hm_w0 = 0;

	//find the relevant breakpoints
	std::vector<double> dconts = rho->getBreakTimes(t0,tt);
	
	if (dconts.size() > 2) {
		//std::cout << "Crossing boundaries!" << std::endl;
	}
	
	//the derivative time integral
	double int_mderiv = 0;
	i = 0;
	for (j = 0; j < dconts.size()-1; j++) {
        if (i >= p->get_length()) { //THIS IS A FUCKING HACK
            break;
        }
		//get the potentials while I'm at it.
		//first the "beginning" potential 
		Hm_w0 += H_wf_r(p->get_traj(i), p->get_time(i), alpha1, alpha2, rho);
		i++;
		//integrate over the interval using the trapezoid rule
		while (p->get_time(i) < dconts[j+1]) {
            if (p->get_traj(i) < 0 || p->get_traj(i) >= PI) {
                return -INFINITY; //Make sure the proposed path is stuck in the right space
            }
			int_mderiv += (dadx_wf_r(p->get_traj(i),p->get_time(i),alpha1,alpha2,rho)+dadx_wf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha2,rho))/2.0
			* (p->get_time(i)-p->get_time(i-1));
			i++;
		}

		//and the last little bit, where I need a left limit
		int_mderiv += (dadx_wf_r(p->get_traj(i),p->get_time(i),alpha1,alpha2,rho,1)+dadx_wf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha2,rho))/2.0
		* (p->get_time(i)-p->get_time(i-1));
		//then the "end" potential
        double tmp = H_wf_r(p->get_traj(i), p->get_time(i), alpha1, alpha2, rho, 1);
        if (tmp != tmp) {
            std::cout << "computation of H_wf_r is nan" << std::endl;
            std::cout << "p->get_length() = " << p->get_length() << " i = " << i << " j = " << j << std::endl;
        }
        Hm_wt += tmp;
	}

	//compute the time integral of the square
	double int_msquare = 0;
	i = 1;
	for (j = 0; j < dconts.size()-1; j++) {
        if (i == p->get_length()) { //THIS IS A FUCKING HACK
            break;
        }
		//integrate over the interval using the trapezoid rule
		while (p->get_time(i) < dconts[j+1]) {
			int_msquare += (a2_wf_r(p->get_traj(i),p->get_time(i),alpha1,alpha2,rho)+a2_wf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha2,rho))/2.0
			* (p->get_time(i)-p->get_time(i-1));
			i++;
		}
        //i--;
        if (i == p->get_length())  {
            std::cout << "i is wrong" << std::endl;
            std::cout << "i = " << i << " p->get_length() = " << p->get_length() << std::endl;
            std::cout << "j = " << j << std::endl;
            std::cout << "dconts = " << std::endl;
            for (int k = 0; k < dconts.size(); k++) {
                std::cout << dconts[k] << " ";
            }
            std::cout << std::endl;
            raise(SIGSEGV);
        }		//and the last little bit, where I need a left limit
		int_msquare += (a2_wf_r(p->get_traj(i),p->get_time(i),alpha1,alpha2,rho,1)+a2_wf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha2,rho))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}
	
	//compute the time integral of the time derivative
	double int_mtime = 0;
	i = 1;
	for (j = 0; j < dconts.size()-1; j++) {
        if (i == p->get_length()) { //THIS IS A FUCKING HACK
            break;
        }
		//integrate over the interval using the trapezoid rule
		while (p->get_time(i) < dconts[j+1]) {
			int_mtime += (dHdt_wf_r(p->get_traj(i),p->get_time(i),alpha1,alpha2,rho)+dHdt_wf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha2,rho))/2.0
			* (p->get_time(i)-p->get_time(i-1));
			i++;
		}
        //and the last little bit, where I need a left limit
		int_mtime += (dHdt_wf_r(p->get_traj(i),p->get_time(i),alpha1,alpha2,rho,1)+dHdt_wf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha2,rho))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}
	
	if (!is_bridge) {
        double gir = (Hm_wt-Hm_w0-1.0/2.0*int_mderiv-1.0/2.0*int_msquare-int_mtime);
        if (gir != gir) {
            std::cout << "ERROR: log_girsanov_wf_r is nan" << std::endl;
            std::cout << "Hm_wt = " << Hm_wt << " Hm_w0 = " << Hm_w0 << " int_mderiv = " << int_mderiv << " int_msquare = " << int_msquare << " int_mtime = " << int_mtime << std::endl;
        }
        return gir;
	} else {
		double gir = (Hm_wt-Hm_w0-1.0/2.0*int_mderiv-1.0/2.0*int_msquare-int_mtime);
		double cond = log_transition_density(x0, xt, taut-tau0);
		return gir + cond;
	}
}

double cbpMeasure::log_girsanov_wfwf(path* p, double alpha1, double alpha2) {
	int i;
	int path_len = p->get_length();
	double x0 = p->get_traj(0);
	double xt = p->get_traj(path_len - 1);
	double t0 = p->get_time(0);
	double tt = p->get_time(path_len -1);
	
	//compute everything for test measure m
	double Hm_wt = H_wfwf(xt,tt,alpha1,alpha2);
	double Hm_w0 = H_wfwf(x0, t0,alpha1,alpha2);
	//compute the time integral of the derivative
	double int_mderiv = 0;
	for (i = 1; i < path_len; i++) {
		//also check if this thing is legal
		if (p->get_traj(i) < 0 || p->get_traj(i) >= PI) {
			return -INFINITY;
		}
		//use trapezoid rule: average the values at t_i, t_{i-1}
		int_mderiv += (dadx_wfwf(p->get_traj(i),p->get_time(i),alpha1,alpha2)+dadx_wfwf(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha2))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}
	//compute the time integral of the square
	double int_msquare = 0;
	for (i = 1; i < path_len; i++) {
		int_msquare += (a2_wfwf(p->get_traj(i),p->get_time(i),alpha1,alpha2)+a2_wfwf(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha2))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}

	return (Hm_wt-Hm_w0-1.0/2.0*int_mderiv-1.0/2.0*int_msquare);
}

double cbpMeasure::log_girsanov_wfwf_r(path* p, double alpha1, double alpha1p, double alpha2, double alpha2p, popsize* rho) {
	int i;
	int j;
	int path_len = p->get_length();

	//get times
	double t0 = p->get_time(0);
	double tt = p->get_time(path_len -1);
	//convert to tau

	
	//compute everything for test measure m
	double Hm_wt = 0;
	double Hm_w0 = 0;
	
	//find the relevant breakpoints
	std::vector<double> dconts = rho->getBreakTimes(t0,tt);
	
	//the derivative time integral
	double int_mderiv = 0;
	i = 0;
	for (j = 0; j < dconts.size()-1; j++) {
        if (i >= p->get_length()) { //THIS IS A FUCKING HACK
            break;
        }
		//get the potentials while I'm at it.
		//first the "beginning" potential 
		Hm_w0 += H_wfwf_r(p->get_traj(i), p->get_time(i), alpha1, alpha1p, alpha2, alpha2p, rho);
		
		i++;
		//integrate over the interval using the trapezoid rule
		while (p->get_time(i) < dconts[j+1]) {
            if (p->get_traj(i) < 0 || p->get_traj(i) > PI) {
                return -INFINITY; //Make sure the proposed path is stuck in the right space
            }
			int_mderiv += (dadx_wfwf_r(p->get_traj(i),p->get_time(i),alpha1, alpha1p,alpha2,alpha2p,rho)+dadx_wfwf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha1p,alpha2,alpha2p,rho))/2.0
			* (p->get_time(i)-p->get_time(i-1));
			i++;
		}
        //and the last little bit, where I need a left limit
		int_mderiv += (dadx_wfwf_r(p->get_traj(i),p->get_time(i),alpha1,alpha1p,alpha2,alpha2p,rho,1)+dadx_wfwf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha1p,alpha2,alpha2p,rho))/2.0
		* (p->get_time(i)-p->get_time(i-1));
		//then the "end" potential
		Hm_wt += H_wfwf_r(p->get_traj(i), p->get_time(i), alpha1,alpha1p, alpha2,alpha2p, rho, 1);
	}
	
	//compute the time integral of the square
	double int_msquare = 0;
	i = 1;
	for (j = 0; j < dconts.size()-1; j++) {
        if (i == p->get_length()) { //THIS IS A FUCKING HACK
            break;
        }
		//integrate over the interval using the trapezoid rule
		while (p->get_time(i) < dconts[j+1]) {
			int_msquare += (a2_wfwf_r(p->get_traj(i),p->get_time(i),alpha1,alpha1p,alpha2,alpha2p,rho)+a2_wfwf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha1p,alpha2,alpha2p,rho))/2.0
			* (p->get_time(i)-p->get_time(i-1));
			i++;
		}
        //and the last little bit, where I need a left limit
        int_msquare += (a2_wfwf_r(p->get_traj(i),p->get_time(i),alpha1,alpha1p,alpha2,alpha2p,rho,1)+a2_wfwf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha1p,alpha2,alpha2p,rho))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}
	
	//compute the time integral of the time derivative
	double int_mtime = 0;
	i = 1;
	for (j = 0; j < dconts.size()-1; j++) {
        if (i == p->get_length()) { //THIS IS A FUCKING HACK
            break;
        }
		//integrate over the interval using the trapezoid rule
		while (p->get_time(i) < dconts[j+1]) {
			int_mtime += (dHdt_wfwf_r(p->get_traj(i),p->get_time(i),alpha1,alpha1p,alpha2,alpha2p,rho)+dHdt_wfwf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha1p,alpha2,alpha2p,rho))/2.0
			* (p->get_time(i)-p->get_time(i-1));
			i++;
		}
        //and the last little bit, where I need a left limit
		int_mtime += (dHdt_wfwf_r(p->get_traj(i),p->get_time(i),alpha1,alpha1p,alpha2,alpha2p,rho,1)+dHdt_wfwf_r(p->get_traj(i-1),p->get_time(i-1),alpha1,alpha1p,alpha2,alpha2p,rho))/2.0
		* (p->get_time(i)-p->get_time(i-1));
	}
	
    double gir = Hm_wt-Hm_w0-1.0/2.0*int_mderiv-1.0/2.0*int_msquare-int_mtime;
    
    if (gir != gir) {
        std::cout << "ERROR: log_girsanov_wfwf_r is nan" << std::endl;
        std::cout << "Hm_wt = " << Hm_wt << " Hm_w0 = " << Hm_w0 << " int_mderiv = " << int_mderiv << " int_msquare = " << int_msquare << " int_mtime = " << int_mtime << std::endl;
    }
	
	return gir;

}


path* wienerMeasure::prop_path(double x0, double t0, double t, std::vector<double>& time_vec) {
	std::vector<double> traj(time_vec.size(),0);
	traj[0] = x0;
	for (int i = 1; i < time_vec.size(); i++) {
		traj[i] = traj[i-1] + random->normalRv(0,sqrt(time_vec[i]-time_vec[i-1]));
	}
	
	path* bm_path = new path(traj, time_vec);
	
	return bm_path;
}

path* wienerMeasure::make_bb_from_bm(path* bm,double u, double v) {
	double T = bm->get_time(bm->get_length()-1);
	double t0 = bm->get_time(0);
	double bT = bm->get_traj(bm->get_length()-1);
	std::vector<double> p;
	for (int i = 0; i < bm->get_length(); i++) {
		double cur_val = (1-(bm->get_time(i)-t0)/(T-t0))*u + (bm->get_time(i)-t0)/(T-t0)*v+bm->get_traj(i)-(bm->get_time(i)-t0)/(T-t0)*bT;
		p.push_back(cur_val);
	}
	path* bb = new path(p, bm->get_time_ref());
	return bb;
}

path* wienerMeasure::prop_bridge(double x0, double xt, double t0, double t, std::vector<double>& time_vec) {
	path* bm;
	path* bb;
	bm = prop_path(0,t0,t,time_vec);
	bb = make_bb_from_bm(bm, x0, xt);
	delete bm;
	return bb;
}

//adapted from Wood (1994), Simulation of the von Mises Fisher distribution and the R package
//IMPORTANT:
//only works for mean vector (0,0,...,0,1), which is okay for this purpose, thank god...
std::vector<double> cbpMeasure::rvMF(double kappa, int d) {
	int i;
	if (kappa == 0) {
		return unifSphere(d);
	} else if (d == 1) {
		double u = random->uniformRv();
		double p = 1.0/(1.0+exp(2.0*kappa));
		std::vector<double> y(1);
		if (u < p) {
			y[0] = -1;
			return y;
		} else {
			y[0] = 1;
			return y;
		}
	} else {
		double W = rW(kappa,d);
		std::vector<double> V = unifSphere(d-1);
		for (i = 0; i < d-1; i++) {
			V[i] = sqrt(1-W*W)*V[i];	
		}
		V.push_back(W);
		return V;
	}
}

//generate a uniform random varaible on the sphere
std::vector<double> cbpMeasure::unifSphere(int d) {
	std::vector<double> y(d,0);
	int i;
	double sum_square = 0;
	for (i = 0; i < d; i++) {
		y[i] = random->normalRv(0,1);
		sum_square += pow(y[i],2);
	}
	double norm = sqrt(sum_square);
	for (i = 0; i < d; i++) {
		y[i] = y[i]/norm;
	}
	return y;
}

//generate W; see Wood (1994)
//CITE THE R FUNCTION
double cbpMeasure::rW(double kappa, int m) {
	double l = kappa;
	double d = m-1;
	double b = (-2.0*l + sqrt(4.0*l*l+d*d))/d;
	double x = (1.0-b)/(1.0+b);
	double c = l*x+d*log(1.0-x*x);
	double u, w, z;
	
	bool done = 0;

	while (!done) {
		z = random->betaRv(d/2.0, d/2.0);
		w = (1.0-(1.0+b)*z)/(1.0-(1.0-b)*z);
		u = random->uniformRv();
		if (l*w+d*log(1.0-x*w)-c >= log(u)) {
			done = 1;
		}
	}
	return w;
}

path* cbpMeasure::prop_bridge(double x0, double xt, double t0, double t, std::vector<double>& time_vec) {
	wienerMeasure myWiener(random);
	int i;
	std::vector<double> u(4,0);
	u[3] = x0;
	double kappa = x0*xt/(t-t0);
	std::vector<double> v = rvMF(kappa,4);
	std::vector<path*> bb_paths(4,NULL);
	for (i = 0; i < 4; i++) {
		bb_paths[i] = myWiener.prop_bridge(u[i], xt*v[i], t0, t, time_vec);
	}
	std::vector<double> b4_traj;
	for (int j = 0; j < bb_paths[0]->get_length(); j++) {
		b4_traj.push_back(0);
		for (i = 0; i < 4; i++) {
			b4_traj[j] += pow(bb_paths[i]->get_traj(j),2);
			if (b4_traj[j] != b4_traj[j]) {
                std::cout << "Failing to propose a BES4 bridge from " << x0 << " to " << xt << " during time interval (" << t0 << ", " << t << ")" << std::endl;
                std::cout << "This likely means that the time vector is getting loopy, possibly due to pileup of points" << std::endl;
                std::cout << "The " << i << "th Brownian bridge between " << u[i] << " and " << xt*v[i] << " is faulty:" << std::endl;
				bb_paths[i]->print_traj(std::cout);
				bb_paths[i]->print_time(std::cout << std::setprecision(20));
			}
		}
		b4_traj[j] = sqrt(b4_traj[j]);
	}
	path* bes4_bridge = new path(b4_traj,time_vec);
	//clean up some memory
	for (i = 0; i < bb_paths.size(); i++) {
		delete bb_paths[i];
	}
	bb_paths.clear();
	return bes4_bridge;
}

//NOTE: parameters are as if in the UNFLIPPED case
path* flippedCbpMeasure::prop_bridge(double x0, double xt, double t0, double t, std::vector<double>& time_vec) {
	cbpMeasure cbp(random);
	path* myPath = cbp.prop_bridge(PI-x0, PI-xt, t0, t,time_vec);
	myPath->flipCbp();
	return myPath;
}

path* wfMeasure::prop_bridge(double x0, double xt, double t0, double t, std::vector<double>& time_vec, double rescale) {
	double dist_from_0 = x0;
	if (xt < x0) dist_from_0 = xt;
	double dist_from_pi = PI-xt;
	if (PI-x0 < PI-xt) dist_from_pi = PI-x0;
	measure* cbp = NULL;
	cbp = new cbpMeasure(random);
//	if (dist_from_0 < dist_from_pi) {
//		cbp = new cbpMeasure(random);
//	} else {
//		cbp = new flippedCbpMeasure(random);
//	}
	if (rescale == -INFINITY) {
		double max_gir = -INFINITY;
		for (int i = 0; i < num_test; i++) {
			path* test_path = cbp->prop_bridge(x0,xt,t0,t,time_vec);
			double gir = cbp->log_girsanov_wf(test_path, 0, 0);
			if (gir > max_gir) {
				max_gir = gir;
			}
			delete test_path;
		}
		rescale = -(log(3) + max_gir);
	}
	bool done = 0;
	path* test_path;
	double gir;
	double accept_prob;
	double u;
    int num_tries = 0;
	while (!done) {
        num_tries++;
		test_path = cbp->prop_bridge(x0, xt, t0, t,time_vec);
		gir = cbp->log_girsanov_wf(test_path, 0, 0);
		accept_prob = rescale + gir;
		if (accept_prob > 0) {
			std::cerr << "Envelope is not sufficient" << std::endl;
			exit(1);
		}
		u = random->uniformRv();
		if (log(u) < accept_prob) {
			done = 1;
		} else {
			delete test_path;
		}
	}
	return test_path;
}

void wfMeasure::invert_path(path* p) {
	for (int i = 0; i < p->get_length(); i++) {
		p->set_traj(inverse_fisher(p->get_traj(i)),i);
	}
}

double cbpMeasure::a(double x, double t) {
	return -1.0/(2.0*x);
}

double cbpMeasure::H(double x, double t) {
	return -log(x)/2.0;
}

double cbpMeasure::dadx(double x, double t) {
	return 1.0/(2.0*x*x);
}


double cbpMeasure::H_wf(double x, double t, double alpha) {
	if (x == 0) {
		return -alpha/2.0;
	} else {
		return log(x)/2.0-1.0/2.0*(alpha*cos(x)+log(sin(x)));
	}
}

double cbpMeasure::a2_wf(double x, double t, double alpha) {
	if (x == 0) {
		return -1.0/6.0*(1+3*alpha);
	} else {
		return -1.0/(4.0*x*x)-1.0/2.0*alpha*cos(x)+1/(4*tan(x)*tan(x))+1.0/4.0*alpha*alpha*sin(x)*sin(x);
	}
}

double cbpMeasure::dadx_wf(double x, double t, double alpha) {
	if (x == 0) {
		return 1.0/6.0*(1+3*alpha);
	} else {
		return -1.0/(2.0*x*x)+1.0/2.0*alpha*cos(x)+1.0/(2.0*sin(x)*sin(x));
	}
}

double cbpMeasure::H_wfwf(double x, double t, double alpha1, double alpha2) {
	return 1.0/2.0*cos(x)*(alpha2-alpha1);
}

double cbpMeasure::a2_wfwf(double x, double t, double alpha1, double alpha2) {
	if (x == 0) {
		return 1.0/2.0*(alpha2-alpha1);
	} else {
		return 1.0/4.0*sin(x)*(alpha1-alpha2)*(sin(x)*(alpha1+alpha2)-2.0/tan(x));
	}
}

double cbpMeasure::dadx_wfwf(double x, double t, double alpha1, double alpha2) {
	return 1.0/2.0*cos(x)*(alpha1-alpha2);
}

double cbpMeasure::H_wf_r(double x, double t, double alpha1, double alpha2, popsize* rho,bool leftLimit) {
	if (x == 0) {
		return -1.0/8.0*rho->getSize(t,leftLimit)*(2*alpha1+alpha2);
	} else {
		return log(x)/2.0-1.0/8.0*(rho->getSize(t,leftLimit)*cos(x)*(2*alpha2+(2*alpha1-alpha2)*cos(x))+4*log(sin(x)));
	}
}

double cbpMeasure::dHdt_wf_r(double x, double t, double alpha1, double alpha2, popsize* rho, bool leftLimit) {
	return -1.0/8.0*rho->getDeriv(t,leftLimit)*cos(x)*(2*alpha2+(2*alpha1-alpha2)*cos(x));
}

double cbpMeasure::a2_wf_r(double x, double t, double alpha1, double alpha2, popsize* rho, bool leftLimit) {
	if (x == 0) {
		return -1.0/6.0*(1.0/rho->getSize(t,leftLimit)+3.0*alpha1);
	} else {
		return 1.0/(16.0*rho->getSize(t,leftLimit))*pow(rho->getSize(t,leftLimit)*sin(x)*(alpha2+(2*alpha1-alpha2)*cos(x))-2*1/tan(x),2)
        - 1.0/(4.0*x*x*rho->getSize(t,leftLimit));
	}
}

double cbpMeasure::dadx_wf_r(double x, double t, double alpha1, double alpha2, popsize* rho, bool leftLimit) {
	if (x == 0) {
		return 1.0/6.0*(1.0/rho->getSize(t, leftLimit)+3.0*alpha1);
	} else {
		return 1.0/2.0*(alpha1*cos(x)+1.0/(sin(x)*sin(x)*rho->getSize(t,leftLimit)))-1.0/(2.0*x*x*rho->getSize(t,leftLimit));
	}
}

double cbpMeasure::H_wfwf_r(double x, double t, double alpha1, double alpha1p, double alpha2, double alpha2p, popsize* rho, bool leftLimit) {
	return 1.0/8.0*rho->getSize(t,leftLimit)*cos(x)*(2*(alpha2-alpha2p)+(2*(alpha1-alpha1p)+alpha2p-alpha2)*cos(x));
}

double cbpMeasure::dHdt_wfwf_r(double x, double t, double alpha1, double alpha1p, double alpha2, double alpha2p, popsize* rho, bool leftLimit) {
	return 1.0/8.0*rho->getDeriv(t,leftLimit)*cos(x)*(2*(alpha2-alpha2p)+(2*(alpha1-alpha1p)+alpha2p-alpha2)*cos(x));
}

double cbpMeasure::a2_wfwf_r(double x, double t, double alpha1, double alpha1p, double alpha2, double alpha2p, popsize* rho, bool leftLimit) {
	return 1.0/64.0*(alpha2p-alpha2+(2*(alpha1p-alpha1)+alpha2-alpha2p)
                     *cos(x))*((-16.0+2.0*rho->getSize(t,leftLimit)*(alpha1+alpha1p)
                     -rho->getSize(t,leftLimit)*(alpha2+alpha2p))
                     *cos(x) + (-2.0*(alpha1+alpha1p)+alpha2+alpha2p)
                     *rho->getSize(t,leftLimit)*cos(3*x)+4*(alpha2+alpha2p)
                     *rho->getSize(t,leftLimit)*sin(x)*sin(x));
}
																																																																	
double cbpMeasure::dadx_wfwf_r(double x, double t, double alpha1, double alpha1p, double alpha2, double alpha2p, popsize* rho, bool leftLimit) {
	return 1.0/2.0*(alpha1p-alpha1)*cos(x);
}


double wfMeasure::a(double x, double t) {
	return 1.0/2.0*(gamma*sin(x)-1/tan(x));
}

double wfMeasure::H(double x, double t) {
	return -1.0/2.0*(gamma*cos(x)+log(sin(x)));
}

double wfMeasure::dadx(double x, double t) {
	return 1.0/2.0*(gamma*cos(x)+1/(sin(x)*sin(x)));
}

double flippedCbpMeasure::a(double x, double t) {
	return 1.0/(2.0*(PI-x));
}

double flippedCbpMeasure::H(double x, double t) {
	return -1.0/2.0*log(2*(PI-x));
}

double flippedCbpMeasure::dadx(double x, double t) {
	return 1.0/(2.0*(PI-x)*(PI-x));
}