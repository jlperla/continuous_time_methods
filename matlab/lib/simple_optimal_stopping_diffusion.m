% Modification of Ben Moll's: http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m
% See notes and equation numbers in 'optimal_stopping.pdf'

% Solves the HJB variational inequality that comes from a general diffusion process with optimal stopping.
%   min{rho v(x) - u(x) - mu(x)v'(x) - sigma(x)^2/2 v''(x), v(x) - S(x)} = 0
%   with a reflecting boundary at a x_min and x_max
%   unless S is very small, and u(x) is very large (i.e. no stopping), the reflecting boundary at x_min is unlikely to enter the solution
%   for a large x_min, it is unlikely to affect the stopping point.

% Does so by using finite differences to discretize into the following complementarity problem:
%   min{rho v - u - A v, v - S} = 0,
%   where A is the discretized intensity matrix that comes from the finite difference scheme and the reflecting barrier at x_min and x_max 

function [results] = simple_optimal_stopping_diffusion(p, settings)

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
        settings.basis_guess = zeros(settings.I,1); %Guess that it never binds?
    end
	%%  Unpack parameters and settings
	rho = p.rho; %Discount rate
	u_x = p.u_x; %utility function
	mu_x = p.mu_x; %Drift function
	sigma_2_x = p.sigma_2_x; %diffusion term sigma(x)^2
	S_x = p.S_x; %payoff function on exit.
	x_min = p.x_min; %Not just a setting as a boundary value occurs here
	x_max = p.x_max; %Not just a setting as a boundary value occurs here.
	
	%Settings for the solution method
	I = settings.I; %number of grid variables for x

	%Create uniform grid and determine step sizes.
	x = linspace(x_min, x_max, I)';

	%% Discretize the operator
	%This is for generic diffusion functions with mu_x = mu(x) and sigma_x = sigma(x)
	mu = mu_x(x); %vector of constant drifts 
	sigma_2 = sigma_2_x(x); %
    
    %Discretize the operator    
    Delta = x(2) - x(1);
    A = discretize_univariate_diffusion(x, mu, sigma_2, false); %Note that this is not checking for absorbing states!
    

	%% Setup and solve the problem as a linear-complementarity problem (LCP)
	%Given the above construction for u, A, and S, we now have the discretized version
	%  min{rho v - u - A v, v - S} = 0,

	%Convert this to the LCP form (see http://www.princeton.edu/~moll/HACTproject/option_simple.pdf)
	%           z >= 0
	%      Bz + q >= 0 
	%   z'(Bz + q) = 0  
	% with the change of variables z = v - S

	u = u_x(x);
	S = S_x(x);
    B = Delta * rho * speye(I) - A; %(6)
    q = -u * Delta + B*S; %(8)

	%% Solve the LCP version of the model 
    %Choose based on the method type.
    if strcmp(settings.method, 'yuval')%Uses Yuval Tassa's Newton-based LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/20952
        %Box bounds, z_L <= z <= z_U.  In this formulation this means 0 <= z_i < infinity
        z_L = zeros(I,1); %(12)
        z_U = inf(I,1);
        settings.error_tolerance = settings.error_tolerance/1000; %Fundamentally different order of magnitude than the others.
        [z, iter, converged] = LCP(B, q, z_L, z_U, settings);
        error = z.*(B*z + q); %(11)
        
     elseif strcmp(settings.method, 'lemke')
        [z,err,iter] = lemke(B, q, settings.basis_guess,settings.error_tolerance, settings.pivot_tolerance);
        error = z.*(B*z + q); %(11)
        converged = (err == 0);

    elseif strcmp(settings.method, 'knitro')
        % Convert the LCP problem to a QCP problem
        % min z' * (B*z - q)
        % s.t. z >= 0
        %      B*z >= -q
        
        % Uses Knitro Tomlab solver
        B = Delta * rho * speye(I) - A; %(6)
        q = -u * Delta + B*S; %(8)
        
        c = zeros(I, 1);
        
        z_iv = zeros(I,1); %initial guess.

        %Box bounds, z_L <= z <= z_U.  In this formulation this means 0 <= z_i < infinity
        z_L = zeros(I,1); %(12)
        z_U = inf(I,1);
        %B*z + q >= 0, b_L <= B*z <= b_U (i.e. -q_i <= (B*z)_i <= infinity)
        b_L = -q;
        b_U = inf(I,1);
        
        Name = 'simple_optimal_stopping_diffusion';
        
        mpec = zeros(I, 6);
        mpec(:, 1) = linspace(1, I, I)';
        mpec(:, 3) = linspace(1, I, I)';
        mpec = sparse(mpec);
        
        
        Prob = lcpAssign(c, z_L, z_U, z_iv, B, b_L, b_U, mpec, Name);
        
        B = Prob.A;
        
        Prob.KNITRO.options.ALG = 3;
        
        Prob.PriLevOpt = 0;   % No output return

        % Solve the LP (with MPEC pairs) using KNITRO:
        Result = tomRun('knitro',Prob,0);
        
        z = Result.x_k(1:I);
        % Original B * original variables, subtract the constants in the
        % constraints to get the result on the c(x) >= 0 form

        dx = B(:,1:I) * z - b_L;
        
        % Should be zero, or very close to zero:
        error = dx(1: I) .* z;
            
        LCP_error = max(abs(error));
        if(LCP_error > settings.error_tolerance)
            if(settings.print_level>0) 
                disp('Failure to converge under Knitro')
            end
            results.success = false;
            return;
        else
            converged = true;
            iter = NaN;
        end
        
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
