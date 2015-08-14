# selection

## Installation

This software is written in C++ and requires the GNU scientific library to run. Download all the files and compile with

```
g++ -O3 -lgsl *.cpp -o sr
```

## Generating allele frequency bridges

To generate a bridge, following the logic of Schraiber et al (2013), use the flag -b, followed by a comma separated list,
`x0,xt,gamma,t`, i.e. the initial frequency, the final frequency, the selection coefficient, and the timespan of the bridge.

For example,

```
./sr -b 0.2,0.6,100,0.2 > myTraj.txt
```

will output a trajectory for an allele with gamma = 100 going from a frequency of 0.2 to 0.6 in 0.2 diffusion time units into the file myTraj.txt. Trajectory files consist of two lines, the first being the allele frequency trajecotry and the second line being the time points.

## Inference from allele frequency time series

Probably you shouldn't do this.

## Other flags that might be relevant

```
-d maximum spacing between time points
-g minimum number of time points
-t number of tests to calibrate rejection sampling algorithm
-e random number seed
```
