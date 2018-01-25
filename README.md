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

1. Sampling times (in years, generations, or units of 2N generations; see below)
2. Sample sizes, in number of (haploid) chromosomes
3. Counts of the derived allele, in number of (haploid) chromsoomes
4. A population size history

It's very important that you are consistent in setting the time scale between items 1 and 4. 

The sample input is a 4-column, white-space-separated file. Each line corresponds to a single sample, which can be just one individual or many individuals from approximately the same time period pooled together. For each sample, each column is

1. The number of derived alleles
2. The sample size (in haploid genomes)
3. The most recent end of the possible age of the sample (i.e. the youngest it could be)
4. The most ancient end of the possible age of the sample (i.e. the oldest it could be)

**NB:** Note that for now sample ages are not sampled during the MCMC, and the sample age is set to the mean of columns 3 and 4

The population size history is a 3-column, white-space-separated file. Each line is one epoch of population size, which can be constant or exponential growth. For each epoch, each column is

1. The population size at the most *ancient* time in the epoch
2. The growth rate of that epoch
3. The most ancient time of the epoch

**VERY IMPORTANT: If sample times are in years, then times in the population size file should be in years. If sample times are in generations, then times in the population size file should be in generations. BE CONSISTENT**

**ALSO IMPORTANT: If times are in units of generations or years, then the growth rates should be in "natural" units. If the times are in diffusion time units, then the growth rates should be scaled by 2N0**

Depending on how you input the data, you may need to specify a generation time and a reference effective population size, N_0, on the command line. If you specified times in units of years you must provide both a generation time and a reference effective population size. If you specified times in units of generations, you must provide a reference effective population size.

On the basis of results from [Jewett et al (2016)](http://biorxiv.org/content/early/2016/04/12/048355.abstract), if no more detailed demography is known, I recommend picking N_0 by computing Waterson's or Tajima's estimators of theta and dividing through by 4 times an estimate of the mutation rate.

The repository includes two sample population size histories, `constant.pop`, which reflects a constant population size, and `horse.all.pop`, which reflects the population history of horses as described by Der Sarkissian et al (2015). If you want to use these in your own analysis, you may need to scale the population sizes. `constant.pop` is scaled for a reference size of N0 = 10000, while `horse.all.pop` is scaled so that the initial population size is 1 (and thus times *must* be specified in diffusion time units).

The essential command line consists of the following flags

```
-D path to data input file
-P path to population size history file
-o output file prefix
-a flag to indicate whether to infer allele age
-G generation time (leave off if times are specified in generations)
-N refernce population size (leave off if times are specified in units of 2N)
```
Note that specifying generation time and population sizes are done with CAPITAL G and CAPITAL N.

This will generate three output files:

1. output_prefix.param: a tab-separated list of samples from the MCMC
2. output_prefix.traj: a list of trajectories sampled from the MCMC
3. output_prefix.time: a list of the times corresponding to each point in each sampled trajectory

Note that the times reported in the output of the MCMC are in units of 2N generations.

You may also wish to control some aspects of the MCMC:

```
-n number of MCMC cycles to run
-f frequency of printing output to the screen
-s frequency of sampling from the posterior
-F fraction of the allele frequency to update during a trajectory update move
```
## Example
For instance, the following command line performs MCMC inference of an allele sampled from a constant population of size 10,000 using the included `constant.pop` file. The dates in `exampleInput.txt` are provided in years, so we must include both a generation time and a base population size: 

```
./sr -D exampleInput.txt -G 25 -N 10000 -n 500000 -d 0.001 -F 20 -f 1000 -s 100 -P constant.pop -e 8067 -a -o output
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
