require(DAAG)
require(expm)

####################MCMC UTILITIES############################

#reads a popsize file
read.popsize = function(popsizefile) {
	params = read.table(popsizefile)
	popSizeFun = function(t) {
		prev_ind = sapply(t,function(s) { min(which(params[,3]<=s)) } )
		prev_time = params[prev_ind,3]
		prev_time[prev_time==-Inf] = 0
		prev_size = params[prev_ind,1]
		rate = params[prev_ind,2]
		popSizes = prev_size*exp(rate*(t-prev_time))
		return(popSizes)
	}
	return(popSizeFun)
}

#make a command string for the MCMC
make_command_string = function(sam_times,sam_sizes,sam_counts,outFile,dt=.001,n=500000,F=20,s=100,f=1000,e=ceiling(10000*runif(1)),popsize_file = "~/Desktop/Selection_Recombination/test.pop",age=TRUE) {
	options(scipen=999)
	command_string = paste("-X",paste(sam_counts,collapse=","),"-N",paste(sam_sizes,collapse=","),"-T",paste(sam_times,collapse=","),"-n",n,"-d",dt,"-F",F,"-f",f,"-s",s,"-P",popsize_file,"-e", e, "-a", "-o",outFile,sep=" ")
	options(scipen=0)
	return(command_string)
}

#reads sampled paths from MCMC
#outname should be the PREFIX, it automatically reads both times and trajectories
read.path = function(outname,n=-1L) {
	traj = readLines(paste(outname,".traj",sep=""),n=n)
	traj = lapply(traj, function(x) {temp = as.numeric(unlist(strsplit(x,split=" "))); temp[2:length(temp)]})
	time = readLines(paste(outname,".time",sep=""),n=n)
	time = lapply(time, function(x) {temp = as.numeric(unlist(strsplit(x,split=" "))); temp[2:length(temp)]})
	return(list(traj=traj,time=time))
}

#plots posterior distribution of paths from MCMC
plot.posterior.paths = function(paths,sam_freqs,sam_times,ylim=c(0,1),truePath=c(),trueTime=c(),dt=.001,plot.ages=T, burnin = 0, xlim=NULL) {
	if (burnin > 0) {
		#ignore paths that are during burnin
		pathsTemp = paths
		paths$time = pathsTemp$time[burnin:length(pathsTemp$time)]
		paths$traj = pathsTemp$traj[burnin:length(pathsTemp$traj)]
		rm(pathsTemp)
	}
	
	#find the oldest age
	oldest = min(sapply(paths$time,min))
	#find the time of the most recent sample 
	present = max(paths$time[[1]])
	#make the vector of times to sample
	post_times = seq(oldest,present,dt)
	#make the matrix to hold the observations from each path
	post_paths = matrix(nrow=length(paths$traj),ncol=length(post_times))
	#loop over every sample path, make it as long as necessary and then run a spline function
	ages = c()
	for (i in 1:length(paths$traj)) {
		ages = c(ages,min(paths$time[[i]]))
		fake_time = seq(oldest-.5,min(paths$time[[i]])-.001,length=1000)
		fake_traj = rep(0,1000)
		cur_traj = c(fake_traj,paths$traj[[i]])
		cur_time = c(fake_time,paths$time[[i]])
		#cur_spline = splinefun(cur_time,cur_traj)
		cur_spline = approxfun(cur_time,cur_traj)
		post_paths[i,] = cur_spline(post_times)
	}	
	path_quantiles = t(apply(post_paths,2,quantile,probs=c(.05,.25,.5,.75,.95)))
	par(mar = c(5, 4, 4, 4) + 0.3)
	five_percent_age = quantile(ages,.05)
	sd_age = sd(ages)
	first_time_ind = min(which(post_times>five_percent_age-sd_age))
	if (is.null(xlim)) {
		matplot(as.numeric(post_times[first_time_ind:length(post_times)]),path_quantiles[first_time_ind:length(post_times),],type="l",lty=1,col=c(3,2,1,2,3),xlab="time",ylab="Allele frequency",ylim=ylim)
	} else {
		matplot(as.numeric(post_times[first_time_ind:length(post_times)]),path_quantiles[first_time_ind:length(post_times),],type="l",lty=1,col=c(3,2,1,2,3),xlab="time",ylab="Allele frequency",ylim=ylim,xlim=xlim)
	}
	points(sam_times,sam_freqs,pch=21,bg="black")
	if (length(truePath)>0 && length(trueTime) > 0) {
		lines(trueTime,truePath,lty=2)
	}
	if (plot.ages) {
		first_nonzero_ind = min(which(sam_freqs>0))
		first_nonzero_time = sam_times[first_nonzero_ind]
		print(c(first_nonzero_ind,first_nonzero_time))
		first_nonzero_post = min(which(post_times>=first_nonzero_time))
		age_dens = density(ages,to=first_nonzero_time)
		age_dens_spline = splinefun(age_dens$x,age_dens$y)
		par(new=T)	
		print("Trying to plot density of times")
		print(c(first_time_ind,first_nonzero_post))
		if (is.null(xlim)) {
			plot(post_times[first_time_ind:length(post_times)],c(age_dens_spline(post_times[first_time_ind:first_nonzero_post]),rep(0,length(post_times)-first_nonzero_post)),type="l",col="blue",axes=FALSE,bty="n",xlab="",ylab="")
		} else {
			plot(post_times[first_time_ind:length(post_times)],c(age_dens_spline(post_times[first_time_ind:first_nonzero_post]),rep(0,length(post_times)-first_nonzero_post)),type="l",col="blue",axes=FALSE,bty="n",xlab="",ylab="",xlim=xlim)
		}
		axis(side=4,at=pretty(range(c(age_dens_spline(post_times[1:first_nonzero_post]),rep(0,length(post_times)-first_nonzero_post)))))
		mtext("Density",side=4,line=3)
	}
	invisible(list(quantiles=path_quantiles,sam_freqs=sam_freqs,sam_times=sam_times,post_times = post_times))
}

##############SIMULATIONS##########################

#Simulate a diploid Wright-Fisher population using an Euler scheme
#popSize is a popsize file read using read.popsize
sim_wf_diploid_popsize = function(a,t_1,t_2,popSize = function(t){1}, alpha2 = 0,h=.5,t_len=1000) {
		wf = vector()
	wf[1] = a
	t = seq(t_1,t_2,length=t_len)
	for (i in 2:length(t)) {
		wf[i] = wf[i-1]+alpha2*wf[i-1]*(1-wf[i-1])*(wf[i-1]+h*(1-2*wf[i-1]))*(t[i]-t[i-1])+sqrt(wf[i-1]*(1-wf[i-1])/popSize(t[i-1]))*sqrt(t[i]-t[i-1])*rnorm(1,0,1)
		wf[i] = min(wf[i],1)
		wf[i] = max(wf[i],0)
		if (wf[i] == 0) {
			wf = c(wf,rep(0,t_len-i))
			break
		} else if (wf[i] == 1) {
			wf = c(wf,rep(1,t_len-i))
			break
		}
	}
	return(rbind(t,wf))
}

#generates sample data from a single path
sample_data_from_path = function(path,sample_times,sample_sizes) {
	#nb: path[1,] is the time, path[2,] is the trajectory
	first_in_t = max(which(sample_times < path[1,1]))+1
	sample_inds = sapply(sample_times[first_in_t:length(sample_times)], function(s){max(which(path[1,]<=s))})
	sample_counts = c(rep(0,first_in_t-1), rbinom(length(sample_inds),prob=path[2,sample_inds],size=sample_sizes))
	freq = c(rep(0,first_in_t-1), path[2,sample_inds])
	return(list(times=sample_times,counts=sample_counts,sizes=sample_sizes,freq=freq))
}

#completely dumb rejection sampler for the age of the allele
#requires ages to come from a finite span of time between ancient and recent
rejection_sample_age = function(n,popSize,ancient,recent=0,M=1) {
	cur_sam = 0
	dat = numeric(n)
	while(cur_sam < n) {
		test = runif(1,ancient,recent)
		u = runif(1)
		if (u<popSize(test)/(M*dunif(test,ancient,recent))) {
			cur_sam = cur_sam + 1
			dat[cur_sam] = test
		}
	}
	return(dat)
	
}

#generates sample data by simulating and sampling alleles
#alpha2, h, and t_1 (which is the age) MUST be vectors of length n
generate_sample_data = function(n,sample_times,sample_sizes, a,t_1,t_2,alpha2,h,t_len=1000, one_nonzero = F, print_i = F, popSize=function(t){1}) {
	#generates n sets of sampling data where samples of size sample_sizes were drawn at sample_times
	#first, simulate the sampling data
	if (is.unsorted(sample_times)) {
		stop("Times are unsorted")
	}
	paths = list()
	sample_counts = matrix(nrow=n,ncol=length(sample_times))
	for (i in 1:n) {
		if (print_i) {
			print(c(i,alpha2[i],h[i],t_1[i]))	
		}
		seg = FALSE
		while (seg == FALSE) {
			test = sim_wf_diploid_popsize(a=a,t_1=t_1[i],t_2=t_2,alpha2=alpha2[i],h=h[i],t_len=t_len, popSize = popSize)
			if (test[2,ncol(test)] < 1 && test[2,ncol(test)] > 0) {	
				if (t_1[i]<sample_times[1]) {
					first_in_t = 1
				} else {
					first_in_t = max(which(sample_times<t_1[i]))+1
				}
				samples_in_path = first_in_t:length(sample_times)
				sample_inds = sapply(sample_times[samples_in_path],function(s){max(which(test[1,]<=s))})
				cur_counts = c(rep(0,length(sample_times)-length(sample_inds)),rbinom(length(sample_inds),prob=test[2,sample_inds],size=sample_sizes[samples_in_path]))
				if (one_nonzero && sum(cur_counts>0) > 0) {
					seg=TRUE
				} else if (!one_nonzero) {
					seg=TRUE
				} 
			}
		}
		paths[[i]] = test
		sample_counts[i,] = cur_counts
	}
	freq = t(t(sample_counts)/sample_sizes)
	return(list(counts=sample_counts,sizes=sample_sizes,times=sample_times,paths=paths,freq=freq))
}
	

#Generate a command string from a list of sims
#... are arguments to be passed to make_command_string
make_command_string_from_sims = function(sim_data, outPrefix, ...) {
	cmd_string = c()
	if (is.matrix(sim_data$counts)) { 
		for (i in 1:nrow(sim_data$counts)) {
			cmd_string = c(cmd_string, make_command_string(sim_data$times,sim_data$sizes,sim_data$counts[i,],paste(outPrefix,i,sep="_"),...))	
		}
	} else {
		cmd_string = c(cmd_string, make_command_string(sim_data$times,sim_data$sizes,sim_data$counts,outPrefix,...))
	}
	return(cmd_string)
}

#bin the data between a and b into num_bin bins
#last_alone = TRUE means the last one is its own bin
bin_data = function(sim_data, a=-.1, b=0, num_bin = 4, bins = NULL, last_alone = TRUE, remove_empty = TRUE) {
	if (!is.vector(bins)) {
		bins = seq(a,b,len=num_bin+1)
	} else {
		num_bin = length(bins)-1
	}
	which.bin = .bincode(sim_data$times,bins)
	new_times = (bins[1:(length(bins)-1)]+bins[2:length(bins)])/2
	num_sim = nrow(sim_data$counts)
	new_counts = matrix(0,nrow=nrow(sim_data$counts),ncol=length(new_times))
	new_sizes = rep(0,length(new_times))
	new_freqs = new_counts
	for (i in 1:length(new_times)) {
		cur_times = (which.bin==i)
		new_counts[,i] = rowSums(matrix(sim_data$counts[,cur_times],nrow=num_sim))
		new_sizes[i] = sum(sim_data$sizes[cur_times])
	}
	if (last_alone) {
		last = length(sim_data$times)
		new_times = c(new_times,sim_data$times[last])
		new_sizes = c(new_sizes,sim_data$sizes[last])
		new_counts = cbind(new_counts,sim_data$counts[,last])
		new_sizes[which.bin[last]] = new_sizes[which.bin[last]] - sim_data$sizes[last]
		new_counts[,which.bin[last]] = new_counts[,which.bin[last]] - sim_data$counts[,last]
	}
	new_freqs = t(t(new_counts)/new_sizes)
	new_data = sim_data
	new_data$counts = new_counts
	new_data$sizes = new_sizes
	new_data$freq = new_freqs
	new_data$times = new_times
	new_data$lower = bins[1:num_bin]
	new_data$upper = bins[2:(num_bin+1)]
	if (last_alone) {
		new_data$lower = c(new_data$lower,sim_data$times[last])
		new_data$upper = c(new_data$upper,sim_data$times[last])
	}
	return(new_data)
}

make_input_matrix_from_sims = function(sim_data,lower=sim_data$times,upper=sim_data$times) {
	inFiles = list()
	empty_bins = sim_data$sizes==0
	num_bins = length(sim_data$sizes[!empty_bins])
	for (i in 1:nrow(sim_data$counts)) {
		curInput = matrix(nrow=num_bins,ncol=4)
		curInput[,1] = sim_data$counts[i,!empty_bins]
		curInput[,2] = sim_data$sizes[!empty_bins]
		curInput[,3] = lower[!empty_bins]
		curInput[,4] = upper[!empty_bins]
		inFiles[[i]] = curInput
	}
	return(inFiles)
}

################HMM LIKELIHOOD#####################

#computes the likelihood
wf_iterate_likelihood_diploid = function(sample_count,sample_size,sample_time,alpha,h,age,N=100) {
	#convert times into generations
	sample_time_gen = floor(2*N*sample_time)
	#convert age to generations
	age_gen = floor(2*N*age)
	#convert alpha to s
	s = alpha/(2*N)
	#make the wf matrix for the parameters
	wf_matrix = construct_wf_matrix_diploid(s,h,N)

	
	#figure out which samples are older than the allele age
	first_sample = min(which(sample_time_gen>age_gen))
	if (first_sample != 1 && sample_count[first_sample-1]>0) {
		print("ERROR: sample with more than 0 copies of the derived allele before age")
		return(0)
	}

	#precompute all necessary transition matrices
	uniqueBetween = unique(c(diff(sample_time_gen),sample_time_gen[first_sample]-age_gen))
	wfTrans = lapply(uniqueBetween,function(s){wf_matrix%^%s})
	
	#now it's just a hidden markov model
	freqs = 0:(2*N)/(2*N)
	initial_states = c(0,1,rep(0,2*N-1))
	time_between = sample_time_gen[first_sample]-age_gen
	transMat = wfTrans[[which(uniqueBetween==time_between)]]
	like = initial_states%*%t(transMat)*dbinom(sample_count[first_sample],sample_size[first_sample],freqs)
	if (first_sample != length(sample_size)) {
		for (i in (first_sample+1):length(sample_size)) {
			time_between = sample_time_gen[i]-sample_time_gen[i-1]
			transMat = wfTrans[[which(uniqueBetween==time_between)]]
			like = like%*%t(transMat)*dbinom(sample_count[i],sample_size[i],freqs)
		}
	}
	if (any(is.na(like))) {
		like = 1e-300
	}
	return(sum(like))
}

#Builds the transition matrix
construct_wf_matrix_diploid = function(s,h,N) {
	starts = seq(0,2*N,1)
	ends = seq(0,2*N,1)
	wf_matrix = sapply(starts,function(i){dbinom(ends,2*N,eta_diploid(i,s,h,N))})
	return(wf_matrix)
}

eta_diploid = function(i,s,h,N) {
	((1+s)*i^2+(1+s*h)*i*(2*N-i))/((1+s)*i^2+2*(1+s*h)*i*(2*N-i)+(2*N-i)^2)
}

#############TRANSFORMATIONS####################

#Fisher's angular transformation
fisher_wf = function(path) {
	acos(1-2*path)
}

#inverts Fisher's angular transformation
inv_wf = function(path) {
	(1-cos(path))/2
}

#############DIFFUSION FUNCTIONS##################

#infinitesimal mean of the transformed Wright-Fisher diffusion 
mu_wf_diploid = function(x,gam,h) {
	1/4*(gam*sin(x)*(1+(2*h-1)*cos(x))-2/tan(x))
}

#potential function of transformed Wright-Fisher diffusion
H_wf_diploid = function(x,gam,h) {
	-1/8*(gam*cos(x)*(2+(2*h-1)*cos(x))+4*log(sin(x)))
}

#derivative of infinitesimal mean of transformed Wright-Fisher diffusion
dmudx_wf_diploid = function(x,gam,h) {
	1/4*(gam*cos(x)+(2*h-1)*gam*cos(2*x)+2*1/sin(x)^2)	
}

#infinitesimal mean of the transformed Wright-Fisher diffusion 
#relative to Bes(0)
mu_wf_diploid_bes0 = function(x,gam,h) {
	res = 1/4*(gam*sin(x)*(1+(2*h-1)*cos(x))-2/tan(x))- -1/2*1/x
	res[x==0] = 0
	return(res)
}

#squared infinitesimal mean of transformed Wright-Fisher diffusion
#relative to Bes(0)
mu_squared_wf_diploid_bes0 = function(x,gam,h) {
	res = (1/4*(gam*sin(x)*(1+(2*h-1)*cos(x))-2/tan(x)))^2- (-1/2*1/x)^2
	res[x==0] = -1/6*(1+3*gam*h)
	return(res)
}

#potential function of transformed Wright-Fisher diffusion
#relative to Bes(0)
H_wf_diploid_bes0 = function(x,gam,h) {
	res = -1/8*(gam*cos(x)*(2+(2*h-1)*cos(x))+4*log(sin(x))) - -1/2*log(x)
	res[x==0] = -1/8*gam*(1+2*h)
	return(res)
}

#derivative of infinitesimal mean of transformed Wright-Fisher diffusion
#relative to Bes(0)
dmudx_wf_diploid_bes0 = function(x,gam,h) {
	res = 1/4*(gam*cos(x)+(2*h-1)*gam*cos(2*x)+2*1/sin(x)^2)-1/(2*x*x)
	res[x==0] = 1/6*(1+3*gam*h)
	return(res)
}

##############GIRSANOV FUNCTIONS##################

#Likelihood of transformed Wright-Fisher path with two different selection coefficients
girsanov_wfwf = function(path,t_vec, alpha1, alpha2, h1, h2) {
	m1 = H_wf_diploid(path[length(path)],alpha1,h1) - 
		H_wf_diploid(path[1],alpha1,h1) -
		1/2*riemann_integral(dmudx_wf_diploid(path,alpha1,h1),t_vec) -
		1/2*riemann_integral(mu_wf_diploid(path,alpha1,h1)^2,t_vec)
	m2 = H_wf_diploid(path[length(path)],alpha2,h2) -
		H_wf_diploid(path[1],alpha2,h2) -
		1/2*riemann_integral(dmudx_wf_diploid(path,alpha2,h2),t_vec) -
		1/2*riemann_integral(mu_wf_diploid(path,alpha2,h2)^2,t_vec)
	return(m1-m2)
}

#Likelihood of Wright-Fisher path relative to Bes(0)
girsanov_wfbes0_limit = function(path,t_vec,alpha,h) {
	m1 = H_wf_diploid_bes0(path[length(path)],alpha,h) - 
		H_wf_diploid_bes0(path[1],alpha,h) - 
		1/2*riemann_integral(dmudx_wf_diploid_bes0(path,alpha,h),t_vec) - 
		1/2*riemann_integral(mu_squared_wf_diploid_bes0(path,alpha,h),t_vec)
	return(m1)
}