%Using the unit testing framework in matlab.  See https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html
%To run tests:
% runtests %would run all of them in the current directory
% runtests('my_test') %runs just the my_test.m file
% runtests('my_test/my_function_test') %runs only `my_function_test function in `my_test'.

function tests = HJBE_discretized_nonuniform_univariate_test
    tests = functiontests(localfunctions);
end

%This is run at the beginning of the test.  Not required.
function setupOnce(testCase)
    addpath('../lib/');
    testCase.TestData.tolerances.test_tol = 1e-9;    
    testCase.TestData.tolerances.lower_test_tol = 1e-8;    %For huge matrices, the inf norm can get a little different.
    testCase.TestData.tolerances.default_csv_precision = '%.10f'; %Should be higher precision than tolerances.test_tol
end

function simple_value_function_test(testCase)
    tolerances = testCase.TestData.tolerances;
    mu_x = @(x) -0.01 * x;
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    u_x = @(x) exp(x);
    rho = 0.05;
    x_min = 1;
    x_max = 2;
    I = 1500;
    x = logspace(log10(x_min),log10(x_max),I)';
    A = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));
    u = u_x(x);
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    [v, success] = simple_HJBE_discretized_univariate(A, x, u, rho);
    %dlmwrite(strcat(mfilename, '_1_v_output.csv'), v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    v_check = dlmread(strcat(mfilename, '_1_v_output.csv'));    
    verifyTrue(testCase,norm(v - v_check, Inf) < tolerances.test_tol, 'v value no longer matches');
    verifyTrue(testCase, success==true, 'unsuccesful');
    
    %Solve with a uniform grid and check if similar after interpolation
    x_2 = linspace(x_min, x_max, 3 * I)'; %Twice as many points to be sure.
    A_2 = discretize_univariate_diffusion(x_2, mu_x(x_2), sigma_2_x(x_2));
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    [v_2, success] = simple_HJBE_discretized_univariate(A_2, x_2, u_x(x_2), rho);

    %Make sure within range.  This seems too large.
    verifyTrue(testCase, norm(interp1(x_2, v_2, x) - v,Inf) < 0.02, 'Not within range of interpolation');    
end

function GBM_adding_points_test(testCase)
% Notice on csv files: 
%1._addpoint is for adding 20 points to the original grid;
%2._22_v_output is result from randomly shifting existing points on the
%grid
%3._base_v_output is result for 1500 points uniform grid;


    tolerances = testCase.TestData.tolerances;
    mu_x = @(x) -0.01 * x;
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    u_x = @(x) exp(x);
    rho = 0.05;
    x_min = 0.1;
    x_max = 3;
    I = 1500;
    I_extra = 20;
    
    %Linspace then merge in extra points
    x_base = linspace(x_min, x_max, I)';
    x_extra = linspace(x_min, x_max, I_extra)';
    x = unique([x_base; x_extra], 'sorted');
    
    %Results from x_base, the uniform grid
    
    A = discretize_univariate_diffusion(x_base, mu_x(x_base), sigma_2_x(x_base));
    u = u_x(x_base);
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    [v_base, success] = simple_HJBE_discretized_univariate(A, x_base, u, rho);
    % This writes the baseline uniform results for v
    %dlmwrite(strcat(mfilename, '_base_v_output.csv'), v_base, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    % Results from x_add
    
    A_a = discretize_univariate_diffusion(x, mu_x(x), sigma_2_x(x));
    u = u_x(x);
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    [v_a, success] = simple_HJBE_discretized_univariate(A_a, x, u, rho);
    % This writes the nonuniform results for v after adding points
    %dlmwrite(strcat(mfilename, '_addpoint_v_output.csv'), v_a, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    
    % Check whether addpoint v is close to uniform v
    v_check = dlmread(strcat(mfilename, '_base_v_output.csv'));    
    verifyTrue(testCase,norm(interp1(x, v_a, x_base) - v_check, Inf) < 1e-3, 'v value no longer matches');
    verifyTrue(testCase, success==true, 'unsuccesful');
    % Check by uniform grids
    
%    x = linspace(x_min, x_max, nn)';    % Data saved in HJBE_discretized_nonuniform_univarite_test_22_v_output.csv

%     % Add only one point to a uniform grid
%     x_base = linspace(x_min, x_max, nn - 1)';
%     index = 3;
%     while length(x_base) ~= nn && index < 19
%         x_base = unique([x_base; x_extra(index)], 'sorted');
%         index = index + 1;    
%     end
%     
%     if length(x_base) == nn
%         x = x_base;
%     else
%         print('Fail to construct a new x.');
%     end
      

end

function NUF_shift_point_test(testCase)
    % Experiment1: Shift most points in from a uniform grid by a tiny number 
    % This should be compared to the results generated from uniform grid
    tolerances = testCase.TestData.tolerances;
    mu_x = @(x) -0.01 * x;
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    u_x = @(x) exp(x);
    rho = 0.05;
    x_min = 0.1;
    x_max = 3;
    I = 1500;

    x_base = linspace(x_min, x_max, I)';
    shifts = (rand(I, 1) - 0.5) / (I * 10e4);
    shifts(1, 1) = 0;
    shifts(end, 1) = 0;
    shifts(601: 610, 1) = zeros(10, 1);
    shifts(1001: 1010, 1) = zeros(10, 1);
    x_s = x_base + shifts;

    A = discretize_univariate_diffusion(x_s, mu_x(x_s), sigma_2_x(x_s));
    u = u_x(x_s);
    
    %Solve with nonuniform grid with random shifts epsilon
    [v, success] = simple_HJBE_discretized_univariate(A, x_s, u, rho);
    %dlmwrite(strcat(mfilename, '_22_v_output.csv'), v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    v_check = dlmread(strcat(mfilename, '_22_v_output.csv'));
    
    %Solve with a uniform grid before using the base.
    %x_2 = linspace(x_min, x_max, 3 * I)'; %Twice as many points to be sure.
    A_2 = discretize_univariate_diffusion(x_base, mu_x(x_base), sigma_2_x(x_base));
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    [v_2, success] = simple_HJBE_discretized_univariate(A_2, x_base, u_x(x_base), rho);
    
    verifyTrue(testCase,norm(interp1(x_base, v_2, x_s) - v_check, Inf) < 1e-6, 'v value no longer matches');
    verifyTrue(testCase, success==true, 'unsuccesful');

end
function NUF_shift_point_2_test(testCase)
    % Experiment2: Shift the nonuniform grid generated by adding point by a
    % tiny number. This should be compared to the results of adding point
    % test
    
    tolerances = testCase.TestData.tolerances;
    mu_x = @(x) -0.01 * x;
    sigma_bar = 0.1;
    sigma_2_x = @(x) (sigma_bar*x).^2;
    u_x = @(x) exp(x);
    rho = 0.05;
    x_min = 0.1;
    x_max = 3;
    I = 1500;
    I_extra = 20;
    
    %Linspace then merge in extra points
    x_base = linspace(x_min, x_max, I)';
    x_extra = linspace(x_min, x_max, I_extra)';
    x = unique([x_base; x_extra], 'sorted');
    
    
    nn = length(x);
    shifts = (rand(nn, 1) - 0.5) / (nn * 10e4);
    shifts(1, 1) = 0;
    shifts(end, 1) = 0;
    shifts(601: 610, 1) = zeros(10, 1);
    shifts(1001: 1010, 1) = zeros(10, 1);
    x_s = x+shifts;
    
    A = discretize_univariate_diffusion(x_s, mu_x(x_s), sigma_2_x(x_s));
    u = u_x(x_s);
    
    %Solve the simple problem: rho v(x) = u(x) + A v(x) for the above process.
    [v_s, success] = simple_HJBE_discretized_univariate(A, x_s, u, rho);
    %dlmwrite(strcat(mfilename, '_2222_v_output.csv'), v, 'precision', tolerances.default_csv_precision); %Uncomment to save again
    v_check = dlmread(strcat(mfilename, '_addpoint_v_output.csv'));    
    verifyTrue(testCase,norm(interp1(x_s, v_s, x)- v_check, Inf) < 1e-6, 'v value no longer matches');
    verifyTrue(testCase, success==true, 'unsuccesful');
    
end