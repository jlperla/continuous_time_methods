% Modification of Ben Moll's: http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m

% Solves the HJB variational inequality that comes from a general diffusion process with optimal stopping.
%   min{rho v(x) - u(x) - mu(x)v'(x) - sigma(x)^2/2 v''(x), v(x) - S(x)} = 0
%   where v'(x_max) = 0, for a given x_max reflecting barrier.
% Does so by using finite differences to discretize into the following complementarity problem:
%   min{rho v - u - A v, v - S} = 0,
%   where A is the discretized intensity matrix that comes from the finite difference scheme and the reflecting barrier at x_max 

function [results] = simple_optimal_stopping_diffusion(p, settings)

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
	Delta = x(2)-x(1); % (1) Would want to use a vector in a non-uniformly spaced grid.
	Delta_2 = Delta^2; %Just squaring the Delta for the second order terms in the finite differences.

	%% Construct sparse A matrix.  TODO: move to a separate file, since general issue.
	%This is for generic diffusion functions with mu_x = mu(x) and sigma_x = sigma(x)
	mu = mu_x(x); %vector of constant drifts 
	sigma_2 = sigma_2_x(x); %

	%Explain the upwinding setup with these X, Y, Z terms.  Given matrices names to make them more clear and map to the code.
    mu_minus = min(mu,0); %General notation of plus/minus.
    mu_plus = max(mu,0); 
    
	X = - mu_minus/Delta + sigma_2/(2*Delta_2);               % (7)
	Y = - mu_plus/Delta + mu_minus/Delta - sigma_2/Delta_2;    % (8)
	Z =  mu_plus/Delta + sigma_2/(2*Delta_2);                % (9)

    %Sparse matrix trick: spdiags takes a vector and an offset.  It returns the vector as sparse a diagonal matrix where the diagonal is offset by the other argument.
    %For example:
    % norm(spdiags([1;2;3], 0, 3, 3)- diag([1 2 3]), Inf) % on the true diagonal, offset 0.  Check that
    % norm(spdiags([2;3;9999], -1, 3, 3)- [0 0 0; 2 0 0; 0 3 0], Inf) %on the diagonal below.  Note that the last element is skipped since only 2 points on off diagonal.
    % norm(spdiags([9999;2;3], 1, 3, 3)- [0 2 0; 0 0 3; 0 0 0], Inf) %on the diagonal above.  Note that the first element is skipped since only 2 points on off diagonal.
    
    %Alternatively this can be done in a single operation to form a tridiagonal matrix by stacking up the arrays, where the 2nd argument is a list of the offsets to apply the columns to)
    %This is equivalent to %A = spdiags(Y, 0, I, I) + spdiags(X(2:I),-1, I, I) + spdiags([0;Z(1:I-1)], 1, I, I);
    A = spdiags([[X(2:I); 0] Y [0; Z(1:I - 1)]], [-1 0 1], I,I);% (10) interior is correct.  Corners will require adjustment    
	
    %Manually adjust the boundary values at the corners.
	A(1,1) = Y(1) + X(1); %Reflecting barrier, (10) and (5)
	A(I,I) = Y(I) + Z(I); %Reflecting barrier,  (10) and (6)

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
	B = rho*speye(I) - A;
	q = -u + B*S; 

	%% Solve the LCP version of the model
	%Uses Yuval Tassa's Newton-based LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/20952
	z_iv = zeros(I,1); %initial guess.

	%Box bounds, z_L <= z <= z_U.  In this formulation this means 0 <= z_i < infinity
	z_L = zeros(I,1);
	z_U = inf(I,1);

	z = LCP(B, q, z_L, z_U, z_iv, settings.display);
	error = z.*(B*z + q);

	LCP_error = max(abs(error));
	if(LCP_error > settings.error_tolerance)
		if(settings.display) 
		disp('Failure to converge')
		end
		results.success = false;
		exit;
	end
	 
	%% Convert from z back to v and plot results
	v = z + S; %calculate value function, unravelling the "z = v - S" change of variables
	
	%% Package Results
    %Discretization results
	results.x = x;
    results.A = A;
    results.S = S;
    
    %Solution
    results.v = v;
	results.success = true;
	results.LCP_error = LCP_error;
end	
