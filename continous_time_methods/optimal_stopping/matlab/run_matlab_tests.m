% Running tests and performance tests.
test_results = runtests('tests','IncludeSubfolders',true);
perf_results = runperf('tests','IncludeSubfolders',true);

fullTable = vertcat(perf_results.Samples);
summaryStats = varfun(@mean,fullTable,...
    'InputVariables','MeasuredTime','GroupingVariables','Name')