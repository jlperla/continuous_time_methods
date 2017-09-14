# Matlab Tests
This folder contains the test-suite for the matlab library and examples in the repository.

## How to Add a Test
* Matlab's unit test framework is based around comments and assertions.  Add as many files in this directory as necessary, as they will all be run by the `run_tests.m`
  * Matlab seems to finds the tests by checking if the filename `test` in it.  So the naming convention of files should end with `_test.m`

* Within a test file then put in an ordered list of the tests with names, e.g. 
```
%% Test 1: MY TEST NAME
% ... code ...
assert(1==1, 'one does not equal one); %Put in checks as an assertion, this would probably pass
assert(0==1, 'zero does not equal one'); %Put in checks as an assertion

%% Test 2: ANOTHER TEST
%Add in other assertions as part of the test.
assert(0==1, 'zero does not equal one'); %Put in checks as an assertion
```

## Tips for Writing Test Conditions
The idea of a test is to comprehensively check that all the results are correct.  Instead of doing this informally during debugging, you should just add the checks into the test file.  At that point, the checks can be run again anytime the code changes.  A few tips:
* When asserting, do not compare floating points with equality, since they may have small differences due to machine precision.  Instead, use a tolerance:
```
tol = 1E-10;
assert(a == a_old ); %Don't use for floating point a and a_old!
assert(abs(a - a_old) < tol); %Do it this way instead.
```
* To compare matrices and vectors, use the `norm`.  The `Inf` norm is the maximum, and is the only one supported by sparse matrices/vectors
```
assert(norm(b - b_old, Inf) < tol); %If b and b_old are vectors 
assert(norm(A - A_old, Inf) < tol); %Also works for matrices, and even sparse matrices.
```

* If the data to check against is too large to write inline in the test, save it as an external file and load it.  For example,
```
% Assume that a large `A` matrix which was previuosly generated and needs to be checked against.
% It could be written to a file (kept in the tests directory) as `test_1_A_output.csv` with:
% dlmwrite('test_1_A_output.csv', A,'precision','%.10f'); %Saves a csv file with 10 digits precision.  Avoid csvread/etc. in matlab since they have limited precision.

%% Test 1: Check A against the old version
%... calculate new A matrix

%Load the old data to compare against with dlmerad, and put into matrix A_old
A_old = dlmread('test_1_A_output.csv');

%Check they are close to the same with the infinity norm.
assert(norm(A - A2, Inf) < tol, 'A has changed compared to the old version')
```
* Save files as CSV where possible, rather than storing matlab files.  The reason is to make it easier to examine changes in `git` and to use the same test file for different languages in a `python`, `julia`, or `C++` port of the algorithm.  To make that sane, use a naming convention for variables.
* Naming convention for output files to compare:
  * GIven a test called `MYTEST.m`, within the `Test 1` section, and with variable `MYVAR, call the file: `MYTEST_1_MYVAR_output.csv`
* For sparse matrices, `dlmwrite` won't work directly, as it only stores dense matrices.  To get around this, you will need to convert the matrix to a sparse format and then convert back when loading. For example, see this roundtrip of saving and storing.
```
 A = 2.0101 * speye(2); %From some sparse matrix
 [indices_i, indices_j, val_ij] = find(A); %Trick to get indicies and values in vectors
 dlmwrite('test_sparse_1_A_output.csv', [indices_i indices_j val_ij],'precision','%.10f'); %Stores these as columns. 
 A_sparse = dlmread('test_sparse_1_A_output.csv');
 A_new = sparse(A_sparse(:,1), A_sparse(:,2), A_sparse(:,3));
 assert(norm(A - A_new, Inf) < test_tol);
 