% Running unit and regression tests and performance tests for all continuous-time methods
addpath('../lib');
clc;
disp('*******Running all Tests*****');
test_results = runtests(pwd,'IncludeSubfolders',true);
