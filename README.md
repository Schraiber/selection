# selection

## Installation

This software is written in C++ and requires the GNU scientific library to run. Download all the files and compile with

```
g++ -O3 -lgsl *.cpp -o sr
```
On some systems such as gcc, the order of arguments matters, in which case
compiling with
```
g++ -O3 *.cpp -lgsl -lgslcblas -lm -o sr
```
may address linker errors.

## Generating allele frequency bridges

To generate a bridge, following the logic of Schraiber et al (2013), use the flag -b, followed by a comma separated list,
`x0,xt,gamma,t`, i.e. the initial frequency, the final frequency, the selection coefficient, and the timespan of the bridge.

For example,

```
./sr -b 0.2,0.6,100,0.2 > myTraj.txt
```

will output a trajectory for an allele with gamma = 100 going from a frequency of 0.2 to 0.6 in 0.2 diffusion time units into the file myTraj.txt. Trajectory files consist of two lines, the first being the allele frequency trajecotry and the second line being the time points.  
A `-R` option causes `sr` to output the trajectory in tsv format.

## Inference from allele frequency time series

The basic input for generating an inference of selection coefficients and allele ages from an allele frequency time series are

1. Sampling times, in units of 2N_0 generations
2. Sample sizes, in number of chromosomes
3. Counts of the derived allele, in number of chromsoomes
4. A population size history

It's very important that you are consistent in setting the time scale between items 1 and 4. On the basis of results from [Jewett et al (2016)](http://biorxiv.org/content/early/2016/04/12/048355.abstract), if no more detailed demography is known, I recommend picking N_0 by computing Waterson's or Tajima's estimators of theta and dividing through by 4 times an estimate of the mutation rate.

The population size history is a 3-column white-space-separated file. Each line is one epoch of population size, which can be constant or exponential growth. For each epoch, each column is

1. The population size at the most *ancient* time in the epoch, relative to N_0
2. The growth rate of that epoch, scaled by 2N_0
3. The most ancient time of the epoch, in units of 2N_0 generations

The repository includes two sample population size histories, `constant.pop`, which reflects a constant population size, and `horse.all.pop`, which reflects the population history of horses as described by Der Sarkissian et al (2015).

The essential command line consists of the following flags

```
-X comma separated list of derived allele counts at each sampling time
-N comma separated list of sample sizes at each sampling time
-T comma separated list of sampling times
-P path to population size history file
-o output file prefix
```

This will generate three output files:

1. output_prefix.param: a tab-separated list of samples from the MCMC
2. output_prefix.traj: a list of trajectories sampled from the MCMC
3. output_prefix.time: a list of the times corresponding to each point in each sampled trajectory

You may also wish to control some aspects of the MCMC:

```
-n number of MCMC cycles to run
-f frequency of printing output to the screen
-s frequency of sampling from the posterior
-F fraction of the allele frequency to update during a trajectory update move
```

For instance, the following command line performs MCMC inference of an allele sampled from a constant population:

```
./sr -X 0,0,2,11,17 -N 20,20,20,20,20 -T -0.6,-0.44994994994995,-0.2997997997998,-0.14964964964965,0 -n 500000 -d 0.001 -F 20 -f 1000 -s 100 -P constant.pop -e 8067 -a -o output
```

## Not inferring allele age

In some cases, it is not as desirable to infer the age of the allele, for example if your allele was ascertained due to being high frequency. In that case, adding inference of allele age will bias your results toward an inference of selection. To deal with this situation, we also have a mode that applies a uniform prior to the population allele frequency at the most ancient sampling time. This can be executed by simply **leaving off** the `-a` flag.

## Other flags that might be relevant

```
-d maximum spacing between time points
-g minimum number of time points
-t number of tests to calibrate rejection sampling algorithm
-e random number seed
```

## Analysis of output

The file `path_utilities.r` has several functions that can be used to visualize output, and also simulate data. Of particularly interest is the function plot.posterior.paths which will make posterior path figures along the lines of Figure 6 in the paper describing the method,
