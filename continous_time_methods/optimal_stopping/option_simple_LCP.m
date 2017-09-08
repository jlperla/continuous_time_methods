% Modification of Ben Moll's: http://www.princeton.edu/~moll/HACTproject/option_simple_LCP.m

% Solves the HJB variational inequality that comes from a general diffusion process with optimal stopping.
%   min{rho v(x) - u(x) - mu(x)v'(x) - sigma(x)^2/2 v''(x), v(x) - S(x)} = 0
%   where v'(x_max) = 0, for a given x_max reflecting barrier.
% Does so by using finite differences to discretize into the following complementarity problem:
%   min{rho v - u - A v, v - S} = 0,
%   where A is the discretized intensity matrix that comes from the finite difference scheme and the reflecting barrier at x_max 

% This version is virtually identical, but has more annotations and comments to aid in learning the method, and is setup to allow specification of more general diffusion processes.
clc; close all;
tic;

%% Parameters and Settings
%Parameters
rho = 0.05; %Discount rate
mu_bar = -0.01; %Drift.  Sign changes the upwind direction.
sigma_bar = 0.01; %Variance
S_bar = 10.0; %the value of stopping
x_min = 0.0; %TODO: Is there something crucial here to check? v(x_min) < S_bar for example for a boundary value?  Or is there even a binding value?
x_max = 1.0; %Reflecting barrier at x_max.  i.e. v'(x_max) = 0 as a boundary value
gamma = 0.5; %u(x) = x^gamma

%Relevant functions for u(x), S(x), mu(x) and sigma(x) for a general diffusion dx_t = mu(x) dt + sigma(x) dW_t, for W_t brownian motion
u_x = @(x) x.^gamma; %u(x) = x^gamma in this example
S_x = @(x) S_bar.*ones(numel(x),1); %S(x) = S_bar in this example
mu_x = @(x) mu_bar * ones(numel(x),1); %i.e. mu(x) = mu_bar
sigma_x = @(x) sigma_bar*x; %i.e. sigma(x) = sigma_bar x
%TODO: In the old version, why does the sigma_bar multiple the grid x, but not the drift?  Is this a mistake in missing out the drift of the brownian motion term?

%Settings for the solution method
N_x = 1000; %number of grid variables for x
display_iterations = false;
error_tolerance = 10^(-6);

%Create uniform grid and determine step sizes.
x = linspace(x_min, x_max, N_x)';
dx = x(2)-x(1);  Would want to use a vector in a non-uniformly spaced grid.
dx_2 = dx^2; %Just squaring the dx for the second order terms in the finite differences.

%% Construct sparse A matrix
%This is for generic diffusion functions with mu_x = mu(x) and sigma_x = sigma(x)
mu = mu_x(x); %vector of constant drifts 
sigma_2 = sigma_x(x).^2; %

%TODO: Explain the upwinding setup with these X, Y, Z terms.  Given matrices names to make them more clear.
X = - min(mu,0)/dx + sigma_2/(2*dx_2);
Y = - max(mu,0)/dx + min(mu,0)/dx - sigma_2/dx_2; %TODO: Why doesn't Y have a /2 in it?
Z =  max(mu,0)/dx + sigma_2/(2*dx_2);

%Creates the core finite difference scheme for the given mu and sigma_2 vectors.
A = spdiags(Y, 0, N_x, N_x) + spdiags(X(2:N_x),-1, N_x, N_x) + spdiags([0;Z(1:N_x-1)], 1, N_x, N_x); %TODO: Better explain the finite difference scheme here with the upwind
A(N_x,N_x)= Y(N_x) + sigma_2(N_x)/(2*dx_2); A(N_x,N_x-1) = X(N_x); %TODO: Better explain how this implments the right hand boundary
%TODO: Explain the left hand boundary value here?  Should we be making sure that v(x_min) <= S, for example...


%TODO: Constructing variations on the creation of A for exposition:
  % Variation 1: construct the A assuming that mu < 0 (i.e., the direction of the finite differences is known a-priori)
  % Variation 2: construct the A with a for loop, essentially adding in each row as an equation.


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
B = rho*speye(N_x) - A;
q = -u + B*S; 

%% Solve the LCP version of the model
%Uses Yuval Tassa's Newton-based LCP solver, download from http://www.mathworks.com/matlabcentral/fileexchange/20952
z_iv = zeros(N_x,1); %initial guess.

%Box bounds, z_L <= z <= z_U.  In this formulation this means 0 <= z_i < infinity
z_L = zeros(N_x,1);
z_U = inf(N_x,1);

z = LCP(B, q, z_L, z_U, z_iv, display_iterations);
error = z.*(B*z + q);

LCP_error = max(abs(error));
if(LCP_error > error_tolerance)
    disp('LCP did not converge')
    exit;
end
 
%% Convert from z back to v and plot results
v = z + S; %calculate value function, unravelling the "z = v - S" change of variables
toc

%Some plots of the results
error = z.*(B*z + q);
p1=plot(x,error);

p2=plot(x,v,x,S,'--','LineWidth',2);
set(gca,'FontSize',16);
legend('v(x)','S(x)','Location','NorthWest');
xlabel('x');
