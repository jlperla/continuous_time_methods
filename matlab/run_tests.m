% Running unit and regression tests and performance tests for all continuous-time methods
addpath('./lib');
clc;
disp('*******Running all Tests*****');
test_results = runtests('tests','IncludeSubfolders',true);
perf_results = runperf('tests','IncludeSubfolders',true);

disp('Summarizing Performance');
fullTable = vertcat(perf_results.Samples);
summaryStats = varfun(@mean,fullTable,...
    'InputVariables','MeasuredTime','GroupingVariables','Name')

%Could test here for speed regressions as required.