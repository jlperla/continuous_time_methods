% Running unit and regression tests and performance tests for all continuous-time methods
clc;
test_results = runtests('tests','IncludeSubfolders',true);
perf_results = runperf('tests','IncludeSubfolders',true);

fullTable = vertcat(perf_results.Samples);
summaryStats = varfun(@mean,fullTable,...
    'InputVariables','MeasuredTime','GroupingVariables','Name')