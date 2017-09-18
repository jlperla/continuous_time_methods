# Matlab Tests
This folder contains the test-suite for the matlab library and examples in the repository.  See  https://www.mathworks.com/help/matlab/matlab_prog/write-simple-test-case-with-functions.html for details on Matlab's unit testing framework (which tries to emulate JUnit, etc.)

## Running the test suite
* To run all the tests, use: `run_tests`
* To run a particular file, example: `runtests('KFE_discretized_univariate_test')`
* To run a particular function in a particular file, example: `runtests('KFE_discretized_univariate_test/small_LLS_vs_eigenvalue_test')`
* To run the performance test, `run_performance_tests`

## How to Add a Test to a File
* For an existing file, just add a new function with `_test` at the end.  The unit testing framework will run it as required.
   * For testing the function, useful to run in isoluation (e.g. `runtests('KFE_discretized_univariate_test/small_LLS_vs_eigenvalue_test'))`
   * To check something, use `verifyTrue(testCase, THE_CONDITION_TO_VERIFY, 'The string if failed...')` or anything in https://www.mathworks.com/help/matlab/matlab_prog/types-of-qualifications.html
* When writing the function, keep in mind that `setupOnce(testCase)` runs at the beginning of any test, and `setup(testCase)` runs before every function separately.  (If those functions exist)

## How to Add a New Test File
* Create a file with `_test.m` as the last string of the name.  Matlab's unit testing will then run as part of the test suite.
* If the function was called `MYTESTFILE_test.m`, then create a function in the file
```
function tests = MYTESTFILE_test
    tests = functiontests(localfunctions);
end
```
* Add in a `setupOnce(testCase)` and/or `setup(testCase)` if required.
* Add in specific functions to test with a descriptive name, and `_test` as the end of the function name.  In general, it is a good idea to split tests up into lots of functions.
