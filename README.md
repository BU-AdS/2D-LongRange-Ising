# Long-Range Ising/Phi4 production code base                    1st February 2019

## 2D phi**4 code

For cross-checking purposes, we have also included a 2D Ising
routine with Metropolis, and Wolff and Swenden Wang cluster updates. 

#### Local

The local routines are straight forward. One simply inputs the desired
\mu^2 and \lambda parameters, and two sets of correlation function
data will be produced. One on the temporal direction, and one purely
spatial. The log file will keep a running tab of all the observables.

#### Non-local

The non-local routines use an exact form of the Wolff algorithm as
described here https://arxiv.org/abs/1401.6805. It is in serial, OMP
parallel, and GPU form. One uses the command line input `sigma` to
adjust the strength of the long range interaction, as detaled in the paper.

We provide two forms of long range coupling: simple power law behaviour
(1/r^sigma) and a radial quantisation inspired (cosh(dt) + cos(dtheta))^sigma.

## Compilation

The code may be compiled for either GPU or CPU architecture. Copy the relevant
Makefile.{ARCH} to Makefile and type 'make'. If you compile for GPU, ensure
the Makefile can see your CUDA lib.

## Using the executable

Once the code is sucessfully compiled, a single executabe will be produced.
One can run the execuatble with the default values defined in the Param 
class (found in `util.h`) or by passing command line arguments. If you use 
the command line to pass arguments, we suggest you make use of a shell 
script to run the executable, as we have shown in the template files.

## Dependencies

For running on the BU cluster, load the following modules

## Contributing

Please create a separate git branch for your contribution, and it will be
peer assesed before merging.

