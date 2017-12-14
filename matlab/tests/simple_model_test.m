
function tests = simple_model_test
    tests = functiontests(localfunctions);
end

%This is run at the beginning of the test.  Not required.
function setupOnce(testCase)
    addpath('../lib/');
end
function simple_v_test(testCase)
%% 1. test on v behavior for time changing u and big t grid
% this test checks when T is large and u is moving sufficiently with t, the
% value functions are smooth
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
    plot(x_base,v(:,1),'-o');hold on
    plot(x_base,v(:,50),'-x'); hold on
    plot(x_base,v(:,N-1),'--','Linewidth',2); hold on
    plot(x_base,v(:,N));
    legend('1st t','50th t','99th t','100th t')
    title('value function plot')

    figure()
    plot(x_grid,diff_1,'-o');hold on
    plot(x_grid,diff_s,'--');hold on
    plot(x_grid,diff_N);
    legend('1st to 2nd','N-2 to N-1','N-1 to N')
    title('difference of v plot')
end

function change_r_test(testCase)
%% 2. test for r change and pi change
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
    % u is time varying

    a = 0.1; % this is defining F(0)=a
    a_h = 5.0; % highest time multiplier value
    %u_tx = @(t,x) exp(x).*((t_max-a)/t_max*t+a); % F(t)=(T-a)/T*t+a
    u_tx = @(t,x) exp(x).*(a_h-abs(t-t_max/2)*(a_h-a)/(t_max/2)); % F(t)=a_h - abs(t-T/2)*(a_h-a_l)/(T/2)
    
    % Uncomment if want to compute v_b, saved in '_3_v_output'

    [t_grid_b, x_grid_b] = meshgrid(t_base,x_base); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
    state_permutations_b = [t_grid_b(:) x_grid_b(:)];

    mu_b = bsxfun(mu_tx, state_permutations_b(:,1), state_permutations_b(:,2)); %applies binary function to these, and remains in the correct stacked order.
    sigma_2_b = bsxfun(sigma_2_tx, state_permutations_b(:,1), state_permutations_b(:,2)); %applies binary function to these, and remains in the correct stacked order.

    %Discretize the operator
    [A_b, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t_base, x_base, mu_b, sigma_2_b);

    u_b = bsxfun(u_tx, state_permutations_b(:,1), state_permutations_b(:,2));
    
    r_center=0.08;
    r_low=0.05;
    for i=1:21
        rho(i) = r_center - abs(i-11)*(r_center-r_low)/10;
        v_b = simple_HJBE_discretized_univariate(A_b, state_permutations_b(:,1), u_b, rho(i)); % state_perm need to be in size N*I
        vv_{i} = reshape(v_b,I,N);
        v_T(:,i)=vv_{i}(:,N); % v at last time node
        v_T1(:,i)=vv_{i}(:,N-1); % v at second to last time node
        v_1(:,i)=vv_{i}(:,1); % v at first time node
        v_2(:,i)=vv_{i}(:,2); % second time node
        v_z1(i,:)=vv_{i}(1,:);% v at z=1th
        v_z500(i,:)=vv_{i}(500,:); % v at z=500th
        v_zI(i,:)=vv_{i}(I,:); % at z last point
    end
    
    % How interest rate change affect v(z=1,t)
    figure()
    plot(rho,v_T(1,:)); hold on
    plot(rho,v_T1(1,:)); hold on
    plot(rho,v_1(1,:),'--');hold on
    plot(rho,v_2(1,:),'--');
    legend('Nth v','N-1th v','1st v','2nd v')
    title('How v(z=1,t) change accross r change')
    ylabel('interest rate')
    
    [rho_t_grid,t_rho_grid] = meshgrid(t_base,rho);
    figure()
    surf(rho_t_grid,t_rho_grid,v_z1)
    title('3D plot for v at z=1st point across r and t')
    ylabel('interest rate')
    xlabel('time')
    
    [point_t_grid,t_rho_grid] = meshgrid(t_base,1:21);
    figure()
    surf(point_t_grid,t_rho_grid,v_z1)
    title('3D plot for v at z=1st point across r points(not value) and t')
    ylabel('rgrid point')
    xlabel('time')
end



