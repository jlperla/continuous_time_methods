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
	x_min = p.x_min; %Not just a setting as a boundar value occurs here
	x_max = p.x_max; %Not just a setting as a boundary value occurs here.
	
	%Settings for the solution method
	I = settings.I; %number of grid variables for x

	%Create uniform grid and determine step sizes.
	x = linspace(x_min, x_max, I)';
	Delta = x(2)-x(1); % Would want to use a vector in a non-uniformly spaced grid.
	Delta_2 = Delta^2; %Just squaring the Delta for the second order terms in the finite differences.

	%% Construct sparse A matrix
	%This is for generic diffusion functions with mu_x = mu(x) and sigma_x = sigma(x)
	mu = mu_x(x); %vector of constant drifts 
	sigma_2 = sigma_2_x(x); %

	%TODO: Explain the upwinding setup with these X, Y, Z terms.  Given matrices names to make them more clear and map to the code.
	X = - min(mu,0)/Delta + sigma_2/(2*Delta_2);               % equation (1)
	Y = - max(mu,0)/Delta + min(mu,0)/Delta - sigma_2/Delta_2;    % equation (2)
	Z =  max(mu,0)/Delta + sigma_2/(2*Delta_2);                % equation (3)

	%Creates the core finite difference scheme for the given mu and sigma_2 vectors.
	A = spdiags(Y, 0, I, I) + spdiags(X(2:I),-1, I, I) + spdiags([0;Z(1:I-1)], 1, I, I);    % equation (4)
	%Term Y is on the diagonal except the last row. Every element below the diagonal is constructed by term X, and ones above the diagonal are contructed by term Z.
	
	A(I,I)= Y(I) + sigma_2(I)/(2*Delta_2); A(I,I-1) = X(I); 
    %The middle rows are represented as equation (5). In particular, the first row is written as equation (6) and the last row is equation (7).
	%The last term on the diagonal becomes "- max(mu,0)/Delta + min(mu,0)/Delta - sigma_2/2*Delta_2", which implements that v(x_max+epsilon)=x_max.
	%TODO: Explain the left hand boundary value here?  Should we be making sure that v(x_min) <= S, for example... Verifying that S(x_min) > 0, for example would be enough?
	%The first row of matrix A implements v(x_min-epsilon)=x_min.


	
	% Variation 1: construct the A assuming that mu < 0 (i.e., the direction of the finite differences is known a-priori)
	X_var1 = - mu/Delta + sigma_2/(2*Delta_2);
	Y_var1 = mu/Delta - sigma_2/Delta_2;
	Z_var1 = sigma_2/(2*Delta_2);
	A_var1 = spdiags(Y_var1, 0, I, I) + spdiags(X(2:I), -1, I, I) + spdiags([0;Z(1:I-1)], 1, I, I);
	A_var1(I,I)= Y(I) + sigma_2(I)/(2*Delta_2);
	% Variation 2: construct the A with a for loop, essentially adding in each row as an equation.  Map to exact formulas in a latex document.
	S = zeros(I+2, I+2);
	for i = 1: I
	  x_i = -mu(i)/Delta + sigma_2(i)/(2*Delta_2);  % equation (8)
	  y_i = mu(i)/Delta - sigma_2(i)/Delta_2;       % equation (9)
	  z_i = sigma_2(i)/(2*Delta_2);              % equation (10)
	  S(i+1, i) = x_i;
	  S(i+1, i+1) = y_i;
	  S(i+1, i+2) = z_i;
	end
	S(I+1, I+1) = mu(I)/Delta - sigma_2(I)/(2*Delta_2);  %% equation (11)
	A_var2 = sparse(S(2: I+1, 2: I+1));
	

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
