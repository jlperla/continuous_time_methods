# Variations on Solutions to Optimal Stopping Problems
See Ben Moll's [Heterogeneous Agent Models in Continuous Time](http://www.princeton.edu/~moll/HACTproject.htm) project.  In particular
* [Stopping Time Problem I: Exercising an Option:](http://www.princeton.edu/~moll/HACTproject/option_simple.pdf): with http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m
* [Stopping Time Problem II: Liquid and Illiquid Assets and Fixed Adjustment Costs:](http://www.princeton.edu/~moll/HACTproject/liquid_illiquid_numerical.pdf): with http://www.princeton.edu/~moll/HACTproject/liquid_illiquid_LCP.m

## Simple Example with a Diffusion
While the application here is for a stopping problem with a diffusion, many of the examples here are simply about the discretizaiton step for the HJBE.
* The file [simple_optimal_stopping_diffusion.m](./simple_optimal_stopping_diffusion.m) is a annotated variation on http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m which uses same linear-complementarity problem solver [LCP.m](./LCP.m) as the HACTproject
* The file [simple_optimal_stopping_diffusion_tomlab.m](./simple_optimal_stopping_diffusion_tomlab.m) is a variation which uses tomlab's linear complementarity solver instead.   While Tomlab's algorithms are overkill for this simple problem, they may become required for more complicated setups with more dimensions or nonlinearity.
* [simple_optimal_stopping_diffusion_nonlinear_tomlab.m](./simple_optimal_stopping_diffusion_tomlab.m): sets up the problem as a (potentially) nonlinear complementarity problem.  While this exactly is clearly linear, a control markov process for the diffusion would lead to nonlinearities, and this serves as a test bed for how to write those types of problems.
* [simple_optimal_stopping_diffusion_non_uniform_grid_tomlab.m](./simple_optimal_stopping_diffusion_non_uniform_grid_tomlab.m): variation with a non-uniform grid.  Important for cases where the nonlinearity occur in known areas of the state space.
