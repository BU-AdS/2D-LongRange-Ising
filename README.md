# BU-AdS 2D-AdS production code base                    17th January 2018

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
of the Ising model.

## Dependencies

In order to compile the eigensolver routines, one must have access a copy
of `EGIEN`. Simply edit `Makefile:Line 6` to be your path to `EIGEN`. We
recommend using version 3.3.4 which can be downloaded from
http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz

## Compilation

Once you have ensured that `EIGEN` is visible to `Makefile`, simply
type `make` to compile the executable `ads_graph`. In addition, we have
included options to `make clean` to clean the ojects created by
`Makefile` and `make tar` to create `ads_graph.tar`.

## Using the executable

Once the code is sucessfully compiled, a single executabe will be produced.
One can run the execuatble with the default values defined in the Param 
class (found in `util.h`) or by passing no command line arguments. Note that 
if you change the default values you must recompile the code. If you use 
the command line to pass arguments, we suggest you make use of a shell 
script to run the executable, as we have shown in the example file 
`example.sh`

## Contributing

The `main` function is very lightweight and the majority of the code is
in the header files. We recommend that additions to this code are written
in the same manner, using the header files rather than separate `.cpp`
files.

