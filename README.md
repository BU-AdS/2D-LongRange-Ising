# BU-AdS 2D-AdS production code base                    31st January 2018

## AdS code

This library is a C++ implementation of a 2D phi**4 theory with long range
interactions mediated by a copy of AdS_{2} at each time slice. It has
the following functionalities for the Laplace operator on the space:

1. CG inverter (propagator)
2. Eigensolver
3. Various hypergeometric and debug utilities

In addition, we have implemented a Monte Carlo based algorithm to perform
importance sampling. The algorithm uses Metropolis steps, with Wolff
steps at a given interval. We also incluse some crude ASCII visualisers
and simple routines to caculate correlation functions and observables
of the Ising model. We also include a set of functions to automatically
compute and store the lattice AdS scaling variables, need to scale the
lattice AdS action against the analytic boundary-boundary propagator.
One loop corrections can also be computed.

## Local and non-local 2D Ising code

For cross-checking purposes, we have also included a 2D Ising
routine, with Wolff cluster updates. The spatial extent is set by
the 'Level' parameter, as in the example.sh file. It is the number of
nodes on the outer level. The number of timeslices is set manually.
One simply inputs the desired \mu^2 and \lambda parameters, and two
sets of correlation function data will be produced. One on the temporal
direction, and one purely spatial.

### Local

The local routines are straight forward.

### Non-local

The non-local routines use and exact form of the Wolff algorithm as
described here https://arxiv.org/abs/1401.6805 At present it is in
a very slow serial form, but improvements will come.

## Dependencies

### EIGEN

In order to compile the eigensolver routines, one must have access a copy
of `EGIEN`. Simply edit the Makefile to be your path to `EIGEN`. We
recommend using version 3.3.4 which can be downloaded from
http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz

### GSL
The lattice tuning routines require linear regression routines, found
in `GSL`. Edit the Makefile to be your path to `GSL`. This build is
using v2.4: https://www.gnu.org/software/gsl/

## Compilation

Once you have ensured that `EIGEN` and 'GSL' are visible to `Makefile`,
simply type `make -j {NUMBER OF PROCESSORS YOU HAVE}` to compile the
executable `adsrun`. In addition, we have included options to
`make clean` to clean the ojects created by `Makefile` and `make tar`
to create `ads_graph.tar`.

## Using the executable

Once the code is sucessfully compiled, a single executabe will be produced.
One can run the execuatble with the default values defined in the Param 
class (found in `util.h`) or by passing command line arguments. Note that 
if you change the default values you must recompile the code. If you use 
the command line to pass arguments, we suggest you make use of a shell 
script to run the executable, as we have shown in the example file 
`example.sh`

## Contributing

Please create a separate git branch for your contribution, and it will be
peer assesed before merging.

