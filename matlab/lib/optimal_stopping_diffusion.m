% Modification of Ben Moll's: http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m
% See notes and equation numbers in 'optimal_stopping.pdf'

% Solves the HJB variational inequality that comes from a general diffusion process with optimal stopping.
%   min{rho v(t,x) - u(t,x) - mu(t,x)D_x v(t,x) - sigma(t,x)^2/2 D_xx v(x) - D_t v(t,x), v(t,x) - S(t,x)} = 0
%   with a reflecting boundary at a x_min and x_max
%   unless S is very small, and u(x) is very large (i.e. no stopping), the reflecting boundary at x_min is unlikely to enter the solution
%   for a large x_min, it is unlikely to affect the stopping point.

% Does so by using finite differences to discretize into the following complementarity problem:
%   min{rho v - u - A v, v - S} = 0,
%   where A is the discretized intensity matrix that comes from the finite difference scheme and the reflecting barrier at x_min and x_max 

function [results] = optimal_stopping_diffusion(p, settings)
    t = p.t;
    x = p.x;
    N = length(t);
    I = length(x);  
    
    %% Default settings
    if ~isfield(settings, 'print_level')
        settings.print_level = 0;
    end
    if ~isfield(settings, 'error_tolerance')
        settings.error_tolerance = 1e-12;
    end
    if ~isfield(settings, 'pivot_tolerance')
        settings.pivot_tolerance = 1e-8;
    end    
    if ~isfield(settings, 'method')
        settings.method = 'yuval'; %Default is the Yuval LCP downloaded from matlabcentral
    end
    if ~isfield(settings, 'basis_guess')
        settings.basis_guess = zeros(I*N,1); %Guess that it never binds?
    end
	%%  Unpack parameters and settings
    [t_grid, x_grid] = meshgrid(t, x); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
    state_permutations = [t_grid(:) x_grid(:)];
    mu = bsxfun(p.mu, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
    sigma_2 = bsxfun(p.sigma_2, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
    S = bsxfun(p.S, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
    u = bsxfun(p.u, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
    
     
    %% Discretize the operator
    %TODO: Should be able to nest the time-varying and stationary ones in this, but can't right now.
    if(N > 1)
        A = discretize_time_varying_univariate_diffusion(t, x, mu, sigma_2);    
    else
        A = discretize_univariate_diffusion(x, mu, sigma_2, false);
    end
    
	%% Setup and solve the problem as a linear-complementarity problem (LCP)
	%Given the above construction for u, A, and S, we now have the discretized version
	%  min{rho v - u - A v, v - S} = 0,

	%Convert this to the LCP form (see http://www.princeton.edu/~moll/HACTproject/option_simple.pdf)
	%           z >= 0
	%      Bz + q >= 0 
	%   z'(Bz + q) = 0  
	% with the change of variables z = v - S

    B = p.rho * speye(I*N) - A; %(6)
    q = -u + B*S; %(8)

	%% Solve the LCP version of the model 
    %Choose based on the method type.
    if strcmp(settings.method, 'yuval')%Uses Yuval Tassa's Newton-based LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/20952
        %Box bounds, z_L <= z <= z_U.  In this formulation this means 0 <= z_i < infinity
        z_L = zeros(I*N,1); %(12)
        z_U = inf(I*N,1);
        settings.error_tolerance = settings.error_tolerance/1000; %Fundamentally different order of magnitude than the others.
        [z, iter, converged] = LCP(B, q, z_L, z_U, settings);
        error = z.*(B*z + q); %(11)
        
     elseif strcmp(settings.method, 'lemke')
        [z,err,iter] = lemke(B, q, settings.basis_guess,settings.error_tolerance, settings.pivot_tolerance);
        error = z.*(B*z + q); %(11)
        converged = (err == 0);

    elseif strcmp(settings.method, 'knitro')
        % Uses Knitro Tomlab as a MPEC solver
        c = zeros(I*N, 1); %i.e. no objective function to minimize.  Only looking for feasibility.
        z_iv = zeros(I*N,1); %initial guess.

        %Box bounds, z_L <= z <= z_U.  In this formulation this means 0 <= z_i < infinity
        z_L = zeros(I*N,1); %(12)
        z_U = inf(I*N,1);
        %B*z + q >= 0, b_L <= B*z <= b_U (i.e. -q_i <= (B*z)_i <= infinity)
        b_L = -q;
        b_U = inf(I*N,1);
        
        %Each row in mpec is a complementarity pair.  Require only 2 non-zeros in each row.
        %In mpec, Columns 1:2 refer to variables, columns 3:4 to linear constraints, and 5:6 to nonlinear constraints:
        % mpec = [   var1,var2 , lin1,lin2 , non1,non2  ; ... ];
        %So a [2 0 3 0 0 0] row would say "x_2 _|_  c_3" for the 3rd linear constrant, and c_3 := A(3,:) x
        num_complementarity_constraints = I*N;
        mpec = sparse(num_complementarity_constraints, 6);
        %The first row is the variable index, and the third is the row of the linear constraint matrix.
        mpec(:, 1) = (1:I*N)';  %So says x_i _|_ c_i for all i.
        mpec(:, 3) = (1:I*N)';
              
        %Creates a LCP
        Prob = lcpAssign(c, z_L, z_U, z_iv, B, b_L, b_U, mpec, 'LCP Problem');
        
        %Add a few settings.  Knitro is the only MPEC solver in TOMLAB
        Prob.PriLevOpt = settings.print_level;
        Prob.KNITRO.options.MAXIT = settings.max_iter;
        if ~isfield(settings, 'knitro_ALG')
            Prob.KNITRO.options.ALG = 3; %Knitro Algorithm.  0 is auto, 3 is SLQP        
        else
            Prob.KNITRO.options.ALG = settings.knitro_ALG;
        end
        Prob.KNITRO.options.BLASOPTION = 0; %Can use blas/mkl... might be more useful for large problems.         
        Prob.KNITRO.options.FEASTOL = settings.error_tolerance; %Feasibility tolerance on linear constraints.
     
        % Solve the LP (with MPEC pairs) using KNITRO:
        Result = tomRun('knitro',Prob);
        z = Result.x_k(1:I*N); %Strips out the slack variables automatically added by the MPEC
  
        error = z.*(B*z + q); %(11)
        converged = (Result.ExitFlag == 0);
        iter = Result.Iter;        
    else
        results = NaN;
        assert(false, 'Unsupported method to solve the LCP');        
    end
	
	%% Package Results
    %% Convert from z back to v
    v = z + S; %(7) calculate value function, unravelling the "z = v - S" change of variables
    
    %Discretization results
	results.x = x;
    results.A = A;
    results.S = S;
    
    %Solution
    results.v = v;
	results.converged = converged;
    results.iterations = iter;
	results.LCP_error = max(abs(error));
    results.LCP_L2_error = norm(error,2);
end	
