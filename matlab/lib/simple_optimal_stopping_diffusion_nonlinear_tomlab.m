% Modification of Ben Moll's: http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m

% TODO: After version with simple_optimal_stopping_diffusion is complete, fork here

%This version of the algorithm will create the setup as a nonlinear complementarity problem.
%While the system here is clearly linear and hence must be a better algorithm, this approach may be useful for extending to complicated non-linear setups, or for jointly solving for prices.