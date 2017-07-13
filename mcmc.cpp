/*
 *  mcmc.cpp
 *  Selection_Recombination
 *
 *  Created by Joshua Schraiber on 5/2/13.
 *  Copyright 2013 UC Berkeley. All rights reserved.
 *
 */

#include "mcmc.h"
#include "MbRandom.h"
#include "settings.h"
#include "path.h"
#include "measure.h"
#include "param.h"

#include<iomanip>
#include<fstream>

mcmc::mcmc(settings& mySettings, MbRandom* r) {
	random = r;
	printFreq = mySettings.get_printFreq();
	sampleFreq = mySettings.get_sampleFreq();
	num_gen = mySettings.get_num_gen();
	minUpdate = mySettings.get_grid();
	if (!mySettings.get_linked()) {
		no_linked_sites(mySettings);
	} else {
		
	}
}

void mcmc::no_linked_sites(settings& mySettings) {
	//open files
	std::ofstream paramFile;
	std::ofstream trajFile;
	std::ofstream timeFile;
	std::string paramName = mySettings.get_baseName() + ".param";
	std::string trajName = mySettings.get_baseName() + ".traj";
	std::string timeName = mySettings.get_baseName() + ".time";
	paramFile.open(paramName.c_str());
	trajFile.open(trajName.c_str());
	timeFile.open(timeName.c_str());
    
    //prepare output file
    if (mySettings.get_infer_age()) {
        paramFile << "gen\tlnL\tpathlnL\talpha1\talpha2\tage\tyt\ttuning.alpha1\ttuning.alpha2\ttuning.age\ttuning.yt" << std::endl;
    } else {
        paramFile << "gen\tlnL\tpathlnL\talpha1\talpha2\ty0\tyt\ttuning.alpha1\ttuning.alpha2\ttuning.y0\ttuning.yt" << std::endl;
    }
	
	//initialize wfMeasure
	wfMeasure* curWF = new wfMeasure(random,0);
	wfMeasure* oldWF = NULL;
	curWF->set_num_test(mySettings.get_num_test());	
	
	//initialize helper measures
	wienerMeasure* myWiener = new wienerMeasure(random);
	measure* myCBP; //this way can decide to flip or not during the mcmc
		
	//propose initial sample path between the sampled points
	wfSamplePath* curPath = new wfSamplePath(mySettings,curWF);
	
	//propose an allele age
	double firstAge;
	if (mySettings.get_infer_age()) {
		firstAge = curPath->get_time(0)-1e-6;//-curWF->expected_age((1-cos(curPath->get_traj(0)))/2.0)/10.0;
		//this essentially assumes that there is no population size change between curPath->get_time(0)-1e-6 and curPath->get_time(0)
		//hopefully this is a reasonable thing to think.
		path* firstPath = new path(curWF->fisher(mySettings.get_fOrigin()), curPath->get_traj(0), firstAge, curPath->get_time(0), curWF, mySettings);
		curPath->set_allele_age(firstAge, firstPath, 0);
		//curPath->print_traj(std::cout);
		//curPath->print_time(std::cout);
	}
	
	param_gamma* alpha1 = new param_gamma(mySettings.get_a1start(),random);
	
	param_gamma* alpha2 = new param_gamma(mySettings.get_a2start(),random);
	
	param_path* curParamPath = new param_path(curPath,alpha1,alpha2,random,mySettings);
	
	//initialize the parameter vector
	std::vector<param*> pars;
	pars.push_back(alpha1);
	pars.push_back(alpha2);
	if (!mySettings.get_infer_age()) {
		pars.push_back(new start_freq(curPath->get_traj(0),random,curParamPath));
	} else {
		pars.push_back(new param_age(firstAge, random, curParamPath, mySettings.get_dt(), mySettings.get_grid()));
	}
	pars.push_back(new end_freq(curPath->get_traj(curPath->get_length()-1), random, curParamPath));
	pars.push_back(curParamPath);
    
	//initialize the proposal ratios
	//probably move this somewhere else
	std::vector<double> propChance(0);
	propChance.push_back(mySettings.get_a1prop()); //update alpha1 1
	propChance.push_back(mySettings.get_a2prop()); //update alpha2 1
	propChance.push_back(mySettings.get_ageprop()); //update start/age 2
	propChance.push_back(mySettings.get_endprop()); //update end 2
	propChance.push_back(mySettings.get_pathprop()); //update path 5
//	propChance.push_back(.1); //gamma -> -gamma .1
//	propChance.push_back(.1); //h -> h->1-h .1
	//store as a cdf
	double sum = 0;
	for (int i = 0; i < propChance.size(); i++) {
		sum += propChance[i];
	}
	double cumsum = 0;
	for (int i = 0; i < propChance.size(); i++) {
		cumsum += propChance[i]/sum;
		propChance[i] = cumsum;
	}
	
	//compute starting lnL
	//curlnL = compute_lnL(curPath, curWF, myWiener);
	curlnL = compute_lnL_sample_only(curPath);
	
	//run mcmc
	for (gen = 0; gen < num_gen; gen++) {
		
		std::string state;
		double propRatio = 0;
		double priorRatio = 0;
		double u = random->uniformRv();
		//propose a parameter change
		for (curProp = 0; curProp < propChance.size(); curProp++) {
			if (u < propChance[curProp]) {
				break;
			}
		}
                
		if (curProp < 5) {
			pars[curProp]->increaseProp();
			propRatio = pars[curProp]->propose();
			priorRatio = pars[curProp]->prior();
        }
			
		oldWF = curWF;
		curWF = new wfMeasure(random,pars[0]->get());
		oldlnL = curlnL;
		curlnL = compute_lnL_sample_only(curPath);

		
		double LLRatio = curlnL-oldlnL;
        if (curlnL != curlnL || oldlnL != oldlnL) {
            std::cout << curlnL << " " << oldlnL << std::endl;
            exit(1);
        }
		if (curProp == 0 || curProp == 1 || curProp == 5 || curProp == 6) {
			//need to compute LL ratio due to the new alpha...
			cbpMeasure quickCBP(random);
			LLRatio += quickCBP.log_girsanov_wfwf_r(curPath, pars[0]->getOld(), pars[0]->get(), pars[1]->getOld(), pars[1]->get(), curPath->get_pop());
		}
		double mh = LLRatio+propRatio+priorRatio;
		u = random->uniformRv();
		
		if (gen % printFreq == 0) {
			std::cout << gen << " " << curProp;
			std::cout << std::setprecision(10) <<  " " << oldlnL << " -> " << curlnL << " " << LLRatio << " " << propRatio << " " << priorRatio << " " << mh << " " << log(u) << " ";
		}
		
		if (log(u) < mh) {
			//accept
			if (curProp < 4) {
				pars[curProp]->increaseAccept();
			}
			curPath->set_update_begin(0);
			curPath->set_old_index(-1);
			delete oldWF;
			state = "Accept";
		} else {
			//reject
			if (curProp < 4) {
				pars[curProp]->reset();
			}
			curPath->reset();
			curPath->set_update_begin(0);
			curPath->set_old_index(-1);
			
			delete curWF;
			curWF = oldWF;
			curlnL = oldlnL;
			state = "Reject";
		}
	
		
		if (gen % printFreq == 0) {
			std::cout << state << std::endl;
//			std::cin.get(); 
		}
	
		
//		if (curProp == 1) {
//			std::cin.get();
//		}
		
		int tuningFreq = num_gen/1000;
		if (tuningFreq < 100) {
			tuningFreq = 100;
		}
		if (tuningFreq > 1000) {
			tuningFreq = 1000;
		}
		if (gen % tuningFreq == 0) {
			for (int i = 0; i < 4; i++) {
				pars[i]->updateTuning();
			}
		}
				
		//std::cout << gen << "\t" << curlnL << "\t" << pars[1]->get() << "\t" << curPath->get_traj(0) << "\t" << curPath->get_traj(curPath->get_length()-1) << std::endl;
        
		if (gen % sampleFreq == 0) {
			cbpMeasure testCBP(random);
			double pathlnL = testCBP.log_girsanov_wf_r(curPath, pars[0]->get(), pars[1]->get(), curPath->get_pop(), 0);
			paramFile << gen << "\t" << curlnL << "\t" << pathlnL << "\t" << pars[0]->get() << "\t" << pars[1]->get() << "\t" << pars[2]->get() << "\t" << curPath->get_traj(curPath->get_length()-1) << "\t" << pars[0]->getTuning() << "\t" << pars[1]->getTuning() << "\t" << pars[2]->getTuning() << "\t" << pars[3]->getTuning() <<  std::endl;
			trajFile << gen << " ";
			curPath->print_traj(trajFile << std::setprecision(20));
			timeFile << gen << " "; 
			curPath->print_time(timeFile << std::setprecision(20));
		}
        //REMOVE THIS
        //std::cin.get();
	}
	
	

}

//for now, this computes the lnL of the WHOLE PATH (wrt Wiener measure) and SAMPLES
//could be optimized to only care about updated portions of path?
double mcmc::compute_lnL(wfSamplePath* p, measure* m, wienerMeasure* wm) {
//	if (curProp == 1) {
//		p->print_traj(std::cout);
//		p->print_time(std::cout);
//	}
	//compute dP/dW
	double gir = wm->log_girsanov(p, m, 0, PI);
	//compute sampling probs
	double sample_prob = 0;
	for (int i = 0; i < p->get_num_samples(); i++) {
//		if (curProp == 1) {
//			std::cout << p->get_sampleFreq(i) << " ";
//			double curSample_prob = p->sampleProb(i);
//			std::cout << curSample_prob << " ";
//		}
		sample_prob += p->sampleProb(i);
	}
//	if (curProp == 1) {
//		std::cout << std::endl;
//		std::cout << ((wfMeasure*)m)->get_gamma() << std::endl;
//		std::cout << gir << " " << sample_prob << std::endl;
//		std::cout << std::endl;
//	}
	//compute the probability of a BM goint the same distance
//	double t_span = p->get_time(p->get_length()-1)-p->get_time(0);
//	double y_span = p->get_traj(p->get_length()-1)-p->get_traj(0);
//	double bm_like = 1.0/2.0*log(2*PI*t_span)-(y_span*y_span)/(2.0*t_span);
	
	if (gir != gir) {
		std::cout << "why the fuck is this shit nan? Generation " << gen << ". Proposal " << curProp << std::endl;
		p->print_traj(std::cout);
		p->print_time(std::cout);
		exit(1);
	}
	return gir + sample_prob;// + bm_like;
}

double mcmc::compute_lnL_sample_only(wfSamplePath* p) {
	double sample_prob = 0;
	for (int i = 0; i < p->get_num_samples(); i++) {
		//		if (curProp == 1) {
		//			std::cout << p->get_sampleFreq(i) << " ";
		//			double curSample_prob = p->sampleProb(i);
		//			std::cout << curSample_prob << " ";
		//		}
		sample_prob += p->sampleProb(i);
	}
	
	return sample_prob;
}

