%Using the unit testing framework in matlab.  See https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
%To run tests:
% runtests %would run all of them in the current directory
% runtests('my_test') %runs just the my_test.m file
% runtests('my_test/my_function_test') %runs only `my_function_test function in `my_test'.

function tests = discretize_time_varying_univariate_diffusion_test
    tests = functiontests(localfunctions);
end

%This is run at the beginning of the test.  Not required.
function setupOnce(testCase)
    addpath('../lib/');
    testCase.TestData.tolerances.test_tol = 1e-9;    
    testCase.TestData.tolerances.default_csv_precision = '%.10f'; %Should be higher precision than tolerances.test_tol
end


function nothing_uniform_test(testCase)
   tolerances = testCase.TestData.tolerances;
    
    mu_tx = @(t, x) -0.1 + t + .1*x;%cludge if constant since bsxfun gets confused otherwise
    sigma_bar = 0.1;
    sigma_2_tx = @(t, x) (sigma_bar*x).^2;
    
    %Want to make sure there are no identical spacing.  To aid in verifying the code.
    x = [0.01 .1 .22 .4 .91 1]';
    t = [0.0 1.1 2.4 5.1 6.9]';
    
   [t_grid, x_grid] = meshgrid(t,x); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
   state_permutations = [t_grid(:) x_grid(:)];
  
   %should check that this does the correct permuations in order and calls the underlying function.
   mu = bsxfun(mu_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
   sigma_2 = bsxfun(sigma_2_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
  
   %Discretize the operator
   [A, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t, x, mu, sigma_2);
   
   %dlmwrite(strcat(mfilename, '_11_A_output.csv'), full(A), 'precision', tolerances.default_csv_precision); %Uncomment to save again
   
   A_check = dlmread(strcat(mfilename, '_11_A_output.csv'));    
    
   verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');    
end

function basic_test(testCase)
    %This will create a time-varying setup, and should show that  
    tolerances = testCase.TestData.tolerances;
    
    mu_tx = @(t, x) -0.1 + t + .1*x;%cludge if constant since bsxfun gets confused otherwise
    sigma_bar = 0.1;
    sigma_2_tx = @(t, x) (sigma_bar*x).^2;
    
    %Grid
    x_min = 0.01;
    x_max = 1;
    I = 5;
	t_min = 0.0;
	t_max = 10.0;
	N = 4;
    x = linspace(x_min, x_max, I)';
	t = linspace(t_min, t_max, N)';
    
   [t_grid, x_grid] = meshgrid(t,x); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
   state_permutations = [t_grid(:) x_grid(:)];
  
   %should check that this does the correct permuations in order and calls the underlying function.
   mu = bsxfun(mu_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
   sigma_2 = bsxfun(sigma_2_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
  
   %Discretize the operator
   [A, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t, x, mu, sigma_2);
   
   %dlmwrite(strcat(mfilename, '_1_A_output.csv'), full(A), 'precision', tolerances.default_csv_precision); %Uncomment to save again
   
   A_check = dlmread(strcat(mfilename, '_1_A_output.csv'));    
    
   verifyTrue(testCase,norm(A - A_check, Inf) < tolerances.test_tol, 'A value no longer matches');
end    

function non_time_varying_test(testCase)
    %This will create a non-time-varying setup, and should show that  
    tolerances = testCase.TestData.tolerances;
    
    % This test use mu(x)=mu*x;sigma(x)^2=sigma^2*x^2 and u(x)=exp(x)
    mu_tx = @(t,x) -0.01 * x+0*t;
    mu_x  = @(x) -0.01 * x;
    sigma_bar = 0.1;
    sigma_2_tx = @(t,x) (sigma_bar*x).^2+0*t;
    sigma_2_x  = @(x) (sigma_bar*x).^2;
    u_tx = @(t,x) exp(x)+t*0;
    u_x  = @(x) exp(x);
    rho = 0.05;
    x_min = 0.1;
    x_max = 3;
    t_min = 0.0;
    t_max = 10.0;
    I = 1500;
    N = 10;
    x = linspace(x_min, x_max, I)';
    t = linspace(t_min,t_max, N)';
    
    % generate permutations of x and t, tacking t first then x
    
    [t_grid, x_grid] = meshgrid(t,x); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
    state_permutations = [t_grid(:) x_grid(:)];

    mu = bsxfun(mu_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
    sigma_2 = bsxfun(sigma_2_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
  
    % Solve for A_n and v_n for non time-varying method, use as check
    %A_n = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));
    %dlmwrite(strcat(mfilename, '_2_An_output.csv'), full(A_n), 'precision', tolerances.default_csv_precision); %Uncomment to save again
    A_check = dlmread(strcat(mfilename, '_2_An_output.csv'));
    u_n=u_x(x);
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the non-time
    %varying sprocess
    %v_n = simple_HJBE_discretized_univariate(A_n, x, u_n, rho);
    %dlmwrite(strcat(mfilename, '_2_vn_output.csv'), v_n, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    v_check = dlmread(strcat(mfilename, '_2_vn_output.csv'));    
    
    % Solve for A and v for time-varying method
    [A, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t, x, mu, sigma_2);
    
    u = bsxfun(u_tx, state_permutations(:,1), state_permutations(:,2)); % sp(:,1) is time sp(:,2) is x
    
    [v,success] = simple_HJBE_discretized_univariate(A, state_permutations(:,1), u, rho); % state_perm need to be in size N*I
    %dlmwrite(strcat(mfilename, '_2_v_output.csv'), v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    v_check_1 = dlmread(strcat(mfilename, '_2_v_output.csv'));
    
    % tests for 1. A is same as A_n; 2. A is same across time; 3. v is same
    % as v_n ; 4. v is same across time
    
    % initial test: make sure the v of time varying result didn't change.
    verifyTrue(testCase,norm(v - v_check_1, Inf) < 1e-9, 'v time-varying result no longer matches');
    
    % this checks the matrix A1 in time varying case same as A in
    % non-time-varying case
    verifyTrue(testCase,norm(A(1:I,1:I) + A(1:I,I+1:2*I) - A_check, Inf) < 1e-5, 'A_1 not match with non-time-varying result'); 
    
    % This checks the matrix A1 in time varying case same as A_n, n
    % randomly drawn
    luck=max(floor(rand*N),2);
    indx=(luck-1)*I+1;
    verifyTrue(testCase,norm(A(1:I,1:I) + A(1:I,I+1:2*I) - A(indx:indx+I-1,indx:indx+I-1) - A(indx:indx+I-1,indx+I:indx+2*I-1) , Inf) < 1e-9, 'A_1 not match with A_n'); 
    %verifyTrue(testCase, success==true, 'unsuccesful');
   
    % This checks the value function of t=0 is same as v in
    % non-time-varying result
    verifyTrue(testCase,norm(v(1:I) - v_check, Inf) < 1e-5, 'v_1 not match with non-time-varying result');
    
    % This checks the v_1 in time varying case same as v_n
    luck=max(floor(rand*N),2);
    indx=(luck-1)*I+1;
    verifyTrue(testCase,norm(v(1:I) - v(indx:indx+I-1), Inf) < 1e-9, 'v_1 not match with v_n');

end    

function time_varying_u_test(testCase)
    %This will create a time-varying setup, and should show that  
    tolerances = testCase.TestData.tolerances;
    
    % mu and sig^2 not time varying
    mu_tx = @(t,x) -0.01 * x+0*t;
    sigma_bar = 0.1;
    sigma_2_tx = @(t,x) (sigma_bar*x).^2+0*t;

    
    %Grid
    rho = 0.05;
    x_min = 0.1;
    x_max = 3;
    I = 1000;
	t_min = 0.0;
	t_max = 1.0;
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
   
   v_b = simple_HJBE_discretized_univariate(A_b, state_permutations_b(:,1), u_b, rho); % state_perm need to be in size N*I
   %dlmwrite(strcat(mfilename, '_3_v_output.csv'), v_b, 'precision', tolerances.default_csv_precision); %Uncomment to save again
   
   [t_grid, x_grid] = meshgrid(t,x); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
   state_permutations = [t_grid(:) x_grid(:)];
  
   mu = bsxfun(mu_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
   sigma_2 = bsxfun(sigma_2_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
  
   %Discretize the operator
   [A, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t, x, mu, sigma_2);
   
   u = bsxfun(u_tx, state_permutations(:,1), state_permutations(:,2));
   
   [v,success] = simple_HJBE_discretized_univariate(A, state_permutations(:,1), u, rho); % state_perm need to be in size N*I
   %dlmwrite(strcat(mfilename, '_33_v_output.csv'), v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
   
   v_check = dlmread(strcat(mfilename, '_3_v_output.csv'));    
   
   
   % test whether baseline result changes
   verifyTrue(testCase,norm(v_b - v_check, Inf) < tolerances.test_tol, 'v_baseline value no longer matches');
   % test whether the add point results after interploration is close to
   % baseline
   v_intp = interp2(t_grid, x_grid, reshape(v,size(x_grid,1),size(x_grid,2)), t_grid_b, x_grid_b); % interpolate v to v_b
   verifyTrue(testCase,norm(reshape(v_intp,size(v_intp,1)*size(v_intp,2),1) - v_check, Inf) < 5e-3, 'v_addpoint value not matches');
   % Notice here: the max diff is 0.0038, not sure whether it's acceptable
   % or not
   
end 

function time_varying_mu_test(testCase)
    %This will create a time-varying setup, and should show that  
    tolerances = testCase.TestData.tolerances;

    
    %Grid
    rho = 0.05;
    x_min = 0.1;
    x_max = 3;
    I = 1000;
	t_min = 0.0;
	t_max = 1.0;
	N = 100;
    x_base = linspace(x_min, x_max, I)';
	t_base = linspace(t_min, t_max, N)';
    a = 0.0; % this is defining F(0)=a
    
    % mu is time varying
    mu_tx = @(t,x) -0.01 * x .*((t_max-a)/t_max*t+a);
    sigma_bar = 0.1;
    sigma_2_tx = @(t,x) (sigma_bar*x).^2+0*t;
    
    I_extra = 15;
    N_extra = 15;
    %Linspace then merge in extra points
    x_extra = linspace(x_min, x_max, I_extra)';
    t_extra = linspace(t_min, t_max, N_extra)';
    %t_extra = t_base;% only change x grid
    %x_extra = x_base; % only change t grid
    x = unique([x_base; x_extra], 'sorted');
    t = unique([t_base; t_extra], 'sorted');
        
    % u is not time varying
    
    u_tx = @(t,x) exp(x); % F(t)=(T-a)/T*t+a

    % Uncomment if want to compute v_b, saved in '_3_v_output'
    
   [t_grid_b, x_grid_b] = meshgrid(t_base,x_base); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
   state_permutations_b = [t_grid_b(:) x_grid_b(:)];
  
   mu_b = bsxfun(mu_tx, state_permutations_b(:,1), state_permutations_b(:,2)); %applies binary function to these, and remains in the correct stacked order.
   sigma_2_b = bsxfun(sigma_2_tx, state_permutations_b(:,1), state_permutations_b(:,2)); %applies binary function to these, and remains in the correct stacked order.
  
   %Discretize the operator
   [A_b, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t_base, x_base, mu_b, sigma_2_b);
   
   u_b = bsxfun(u_tx, state_permutations_b(:,1), state_permutations_b(:,2));
   
   v_b = simple_HJBE_discretized_univariate(A_b, state_permutations_b(:,1), u_b, rho); % state_perm need to be in size N*I
   %dlmwrite(strcat(mfilename, '_4_v_output.csv'), v_b, 'precision', tolerances.default_csv_precision); %Uncomment to save again
   
   [t_grid, x_grid] = meshgrid(t,x); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
   state_permutations = [t_grid(:) x_grid(:)];
  
   mu = bsxfun(mu_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
   sigma_2 = bsxfun(sigma_2_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
  
   %Discretize the operator
   [A, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t, x, mu, sigma_2);
   
   u = bsxfun(u_tx, state_permutations(:,1), state_permutations(:,2));
   
   [v,success] = simple_HJBE_discretized_univariate(A, state_permutations(:,1), u, rho); % state_perm need to be in size N*I
   %dlmwrite(strcat(mfilename, '_44_v_output.csv'), v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
   
   v_check = dlmread(strcat(mfilename, '_4_v_output.csv'));    
   
   
   % test whether baseline result changes
   verifyTrue(testCase,norm(v_b - v_check, Inf) < tolerances.test_tol, 'v_baseline value no longer matches');
   % test whether the add point results after interploration is close to
   % baseline
   v_intp = interp2(t_grid, x_grid, reshape(v,size(x_grid,1),size(x_grid,2)), t_grid_b, x_grid_b); % interpolate v to v_b
   verifyTrue(testCase,norm(reshape(v_intp,size(v_intp,1)*size(v_intp,2),1) - v_check, Inf) < 1e-3, 'v_addpoint value not matches');
   % Notice here: the max diff is 0.0038, not sure whether it's acceptable
   % or not
   
end 

function time_varying_both_shift_test(testCase)
    %This will create a time-varying setup, and should show that  
    tolerances = testCase.TestData.tolerances;

    
    %Grid
    rho = 0.05;
    x_min = 0.1;
    x_max = 3;
    I = 1000;
	t_min = 0.0;
	t_max = 1.0;
	N = 200;
    x_base = linspace(x_min, x_max, I)';
	t_base = linspace(t_min, t_max, N)';
    a = 0.0; % this is defining F(0)=a
    
    % mu is time varying
    mu_tx = @(t,x) -0.01 * x .*((t_max-a)/t_max*t+a);
    sigma_bar = 0.1;
    sigma_2_tx = @(t,x) (sigma_bar*x).^2+0*t;
    
    %Some shifts in x and t spaces
    
    x_shift_1 = floor(0.3 * I);
    x_shift_2 = floor(0.7 * I);
    t_shift_1 = floor(0.3 * N);
    t_shift_2 = floor(0.7 * N);

     
    %Linspace then merge in extra points
    nn = length(x_base);
    shifts = (rand(nn, 1) - 0.5) / (nn * 10e4);
    shifts(1, 1) = 0;
    shifts(end, 1) = 0;
    shifts(x_shift_1 + 1: x_shift_1 + 10, 1) = zeros(10, 1);
    shifts(x_shift_2 + 1: x_shift_2 + 10, 1) = zeros(10, 1);
    x = x_base+shifts;
    
    tt = length(t_base);
    shifts = (rand(tt, 1) - 0.5) / (tt * 10e2);
    shifts(1, 1) = 0;
    shifts(end, 1) = 0;
    shifts(t_shift_1 + 1: t_shift_1 + 10, 1) = zeros(10, 1);
    shifts(t_shift_2 + 1: t_shift_2 + 10, 1) = zeros(10, 1);    
    t = t_base+shifts;   
        
    % u is also time varying
    
    u_tx = @(t,x) exp(x).*((t_max-a)/t_max*t+a); % F(t)=(T-a)/T*t+a

    % Uncomment if want to compute v_b, saved in '_3_v_output'
    
   [t_grid_b, x_grid_b] = meshgrid(t_base,x_base); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
   state_permutations_b = [t_grid_b(:) x_grid_b(:)];
  
   mu_b = bsxfun(mu_tx, state_permutations_b(:,1), state_permutations_b(:,2)); %applies binary function to these, and remains in the correct stacked order.
   sigma_2_b = bsxfun(sigma_2_tx, state_permutations_b(:,1), state_permutations_b(:,2)); %applies binary function to these, and remains in the correct stacked order.
  
   %Discretize the operator
   [A_b, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t_base, x_base, mu_b, sigma_2_b);
   
   u_b = bsxfun(u_tx, state_permutations_b(:,1), state_permutations_b(:,2));
   
   v_b = simple_HJBE_discretized_univariate(A_b, state_permutations_b(:,1), u_b, rho); % state_perm need to be in size N*I
   %dlmwrite(strcat(mfilename, '_5_v_output.csv'), v_b, 'precision', tolerances.default_csv_precision); %Uncomment to save again
   
   [t_grid, x_grid] = meshgrid(t,x); %Generates permutations (stacked by t first, as we want) could look at: [t_grid(:) x_grid(:)]
   state_permutations = [t_grid(:) x_grid(:)];
  
   mu = bsxfun(mu_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
   sigma_2 = bsxfun(sigma_2_tx, state_permutations(:,1), state_permutations(:,2)); %applies binary function to these, and remains in the correct stacked order.
  
   %Discretize the operator
   [A, Delta_p, Delta_m, h_p, h_m] = discretize_time_varying_univariate_diffusion(t, x, mu, sigma_2);
   
   u = bsxfun(u_tx, state_permutations(:,1), state_permutations(:,2));
   
   [v,success] = simple_HJBE_discretized_univariate(A, state_permutations(:,1), u, rho); % state_perm need to be in size N*I
   
   v_check = dlmread(strcat(mfilename, '_5_v_output.csv'));    
   
   
   % test whether baseline result changes
   verifyTrue(testCase,norm(v_b - v_check, Inf) < tolerances.test_tol, 'v_baseline value no longer matches');
   % test whether the add point results after interploration is close to
   % baseline
   v_intp = interp2(t_grid, x_grid, reshape(v,size(x_grid,1),size(x_grid,2)), t_grid_b, x_grid_b); % interpolate v to v_b
   verifyTrue(testCase,norm(reshape(v_intp,size(v_intp,1)*size(v_intp,2),1) - v_check, Inf) < 1e-6, 'v_addpoint value not matches');
   % Notice here: the max diff is 1e-7, not sure whether it's acceptable
   % or not
   
end
