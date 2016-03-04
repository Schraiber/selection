#include <iostream>
#include "matrix.h"
#include "MbRandom.h"
#include "measure.h"
#include "path.h"
#include "settings.h"
#include "mcmc.h"
#include "popsize.h"

int main (int argc, char * const argv[]) {
	
	settings mySettings(argc, argv);
	
	MbRandom* r = new MbRandom(mySettings.get_seed());
	
	if (mySettings.get_p()) {
		mySettings.print();
	}
	
	if (mySettings.get_bridge()) {
		wfMeasure myWF(r,0);
		std::vector<double> pars = mySettings.parse_bridge_pars();
        std::cerr << "Creating bridge with (x0, xt, gamma, t) = (" << pars[0] << ", " << pars[1]
            << ", " << pars[2] << ", " << pars[3] << ")" << std::endl;
		myWF.set_num_test(mySettings.get_num_test());
		myWF.set_gamma(pars[2]);
		path* myPath = new path(myWF.fisher(pars[0]),myWF.fisher(pars[1]),0,pars[3],&myWF,mySettings);
		myWF.invert_path(myPath);
		if (mySettings.get_output_tsv()) {
			myPath->print_tsv(std::cout) ;
		} else {
			myPath->print(std::cout);
		}
		delete myPath;
	} else if (mySettings.get_mcmc()) {
		if (mySettings.get_linked()) {
			
		} else {
			mcmc myMCMC(mySettings,r);
		}
	} else {
		std::cout << "No task specified" << std::endl;
		std::cout << "-b x0,xt,gamma,t for bridge path" << std::endl;
		std::cout << "-X x0,x1,...,xn -N n0,n1,...nn -T t0,t1,...tn for mcmc with just allele frequencies" << std::endl;
	}
		
	delete r;
}
