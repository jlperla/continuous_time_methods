# Variations on Solutions to Optimal Stopping Problems
See Ben Moll's [Heterogeneous Agent Models in Continuous Time](http://www.princeton.edu/~moll/HACTproject.htm) project.  In particular
* [Stopping Time Problem I: Exercising an Option:](http://www.princeton.edu/~moll/HACTproject/option_simple.pdf): with http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m
* [Stopping Time Problem II: Liquid and Illiquid Assets and Fixed Adjustment Costs:](http://www.princeton.edu/~moll/HACTproject/liquid_illiquid_numerical.pdf): with http://www.princeton.edu/~moll/HACTproject/liquid_illiquid_LCP.m

## Simple Example with a Diffusion
* The file [option_simple_LCP.m](./option_simple_LCP.m) is a annotated variation on http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m which uses the sparse(?) linear-complementarity problem solver [LCP.m](./LCP.m)
* The file [option_simple_LCP_tomlab.m](./option_simple_LCP_tomlab.m) is a variation which uses tomlab's linear complementarity solver instead.   While Tomlab's algorithms are overkill for this simple problem, they may be required for more complicated ones.
