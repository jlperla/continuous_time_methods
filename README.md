# Examples of Methods for Continuous Time Problems

## Detailed Notes on Methods
* See [operator_discretization_finite_differences.pdf](operator_discretization_finite_differences.pdf) for details on the discretization of the stochastic proceses , and the solution to simple value function and stationary distribution calculuations.
* See [optimal_stopping.pdf](optimal_stopping.pdf) for details on solving optimal stopping problem (given the discretization).


## Structure of folders
* The source for the documentation of all algorithms is in [docs/](docs/). 
* For the [Matlab](matlab/README.md) specific implementation of the algorithms
    * See library functions to implement the methods is in [matlab/lib/](matlab/lib/)
    * See examples which use the library functions in [matlab/examples/](matlab/examples/)
    * All unit and regression tests of the matlab library are in [matlab/tests/](matlab/tests/).  You can use `run_tests.m` in the matlab folder to run all regression tests after making coding changes.

	
## Variations on Solutions to Optimal Stopping Problems
See Ben Moll's [Heterogeneous Agent Models in Continuous Time](http://www.princeton.edu/~moll/HACTproject.htm) project.  In particular
* [Stopping Time Problem I: Exercising an Option:](http://www.princeton.edu/~moll/HACTproject/option_simple.pdf): with http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m
* [Stopping Time Problem II: Liquid and Illiquid Assets and Fixed Adjustment Costs:](http://www.princeton.edu/~moll/HACTproject/liquid_illiquid_numerical.pdf): with http://www.princeton.edu/~moll/HACTproject/liquid_illiquid_LCP.m
