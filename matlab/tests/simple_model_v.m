clear all
addpath('../lib/');

% mu and sig^2 not time varying
mu_tx = @(t,x) -0.01 * x+0*t;
sigma_bar = 0.1;
sigma_2_tx = @(t,x) (sigma_bar*x).^2+0*t;


%Grid
x_min = 0.1;
x_max = 8;
I = 1000;
t_min = 0.0;
t_max = 10.0;
N = 100;
x_base = linspace(x_min, x_max, I)';
t_base = linspace(t_min, t_max, N)';

I_extra = 15;
N_extra = 15;
%Linspace then merge in extra points
x_extra = linspace(x_min, x_max, I_extra)';
t_extra = linspace(t_min, t_max, N_extra)';
%t_extra = t_base;% only change x grid
%x_extra = x_base; % only change t grid
x = unique([x_base; x_extra], 'sorted');
t = unique([t_base; t_extra], 'sorted');

% u is time varying

a = 0.0; % this is defining F(0)=a
u_tx = @(t,x) exp(x).*((t_max-a)/t_max*t+a); % F(t)=(T-a)/T*t+a

% Uncomment if want to compute v_b, saved in '_3_v_output'

[t_grid_b, x_grid_b] = meshgrid(t_base,x_base); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
state_permutations_b = [t_grid_b(:) x_grid_b(:)];

mu_b = bsxfun(mu_tx, state_permutations_b(:,1), state_permutations_b(:,2)); %applies binary function to these, and remains in the correct stacked order.
sigma_2_b = bsxfun(sigma_2_tx, state_permutations_b(:,1), state_permutations_b(:,2)); %applies binary function to these, and remains in the correct stacked order.

%Discretize the operator
[A_b, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t_base, x_base, mu_b, sigma_2_b);

u_b = bsxfun(u_tx, state_permutations_b(:,1), state_permutations_b(:,2));

rho = 0.09;
v_b = simple_HJBE_discretized_univariate(A_b, state_permutations_b(:,1), u_b, rho); % state_perm need to be in size N*I

v = reshape(v_b,I,N);
diff_s = v(:,N-2) - v(:,N-1);
diff_N = v(:,N-1) - v(:,N);  
diff_1 = v(:,1) - v(:,2);

figure()
plot(x_grid,v(:,1),'-o');hold on
plot(x_grid,v(:,50),'-x'); hold on
plot(x_grid,v(:,N-1),'--'); hold on
plot(x_grid,v(:,N));
legend('1st t','50th t','99th t','100th t')
title('value function plot')

figure()
plot(x_grid,diff_1,'-o');hold on
plot(x_grid,diff_s,'--');hold on
plot(x_grid,diff_N);
legend('1st to 2nd','N-2 to N-1','N-1 to N')
title('difference of v plot')