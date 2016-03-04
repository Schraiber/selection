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

## Inference from allele frequency time series

The basic input for generating an inference of selection coefficients and allele ages from an allele frequency time series are

1. Sampling times, in units of 2N_0 generations
2. Sample sizes, in number of chromosomes
3. Counts of the derived allele, in number of chromsoomes
4. A population size history

It's very important that you are consistent in setting the time scale between items 1 and 4. 

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

## Other flags that might be relevant

```
-d maximum spacing between time points
-g minimum number of time points
-t number of tests to calibrate rejection sampling algorithm
-e random number seed
```
