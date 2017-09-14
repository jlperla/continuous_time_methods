% Running unit and regression tests and performance tests for all continuous-time methods
addpath('../lib');
clc;

disp('*******Running Performance Check *****');
perf_results = runperf(pwd,'IncludeSubfolders',true);

if(numel(perf_results) > 0)
    disp('Summarizing Performance');
    fullTable = vertcat(perf_results.Samples);
    summaryStats = varfun(@mean,fullTable,...
        'InputVariables','MeasuredTime','GroupingVariables','Name')
end
%Could test here for speed regressions as required.