# Variations on Solutions to Optimal Stopping Problems
See Ben Moll's [Heterogeneous Agent Models in Continuous Time](http://www.princeton.edu/~moll/HACTproject.htm) project.  In particular
* [Stopping Time Problem I: Exercising an Option:](http://www.princeton.edu/~moll/HACTproject/option_simple.pdf): with http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m
* [Stopping Time Problem II: Liquid and Illiquid Assets and Fixed Adjustment Costs:](http://www.princeton.edu/~moll/HACTproject/liquid_illiquid_numerical.pdf): with http://www.princeton.edu/~moll/HACTproject/liquid_illiquid_LCP.m

## Simple Example with Optimal Stopping of a Diffusion Process
* The file [simple_optimal_stopping_diffusion.m](./simple_optimal_stopping_diffusion.m) is a annotated variation on http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m which uses same linear-complementarity problem solver [LCP.m](./LCP.m) as the HACTproject
* The file [simple_optimal_stopping_diffusion_tomlab.m](./simple_optimal_stopping_diffusion_tomlab.m) is a variation which uses tomlab's linear complementarity solver instead.   While Tomlab's algorithms are overkill for this simple problem, they may become required for more complicated setups with more dimensions or nonlinearity.
