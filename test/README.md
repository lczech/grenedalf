# Test

This directory contains test cases for the estimators implemented in grenedalf, in particuar Theta Pi, Theta Watterson, Tajima's D, and FST. The tests span a wide range of parameters, i.e. different pool sizes, read depths, allele frequencies, and window sizes. However, this is a minimalistic test with clean data, i.e., we assume sufficient read depth and no missing data in the test.

The tests are used for internal testing, but can also be used as a test dataset for other tools as well. The created test data can for instance be useful as a standardized test for evaluating other tools. See the `test_results.csv` table here for the full range of parameters currently tested, as well as the results obtained with grenedalf and PoPoolation.

## Test setup

The creation of test cases and files is implemented in `execute_tests.py`. Then, we run grenedalf as well as our minimalistic independent implementation of the equations in Python, and compare the two to each other. As they both yield the same results, we gain some confidence that the equations are correctly implemented. This is used as an automated test of grenedalf, as we run these tests in our [CI Tests](https://github.com/lczech/grenedalf/actions) in GitHub Actions.

If you want to use the test cases for your own tool tests, have a look at the `test_results.csv` table, as well as the produced `mpileup` and `sync` directories, which contain simple files representing the test cases. Note that these assume error-free data. Hence, if you encounter any issues with grenedalf on more realistic data, please report the [issue](https://github.com/lczech/grenedalf/issues).

The test scripts are also published in the supporting grenedalf manuscript repository, see [here](https://github.com/lczech/grenedalf-paper/tree/master/eval-independent-test). That repository contains the scripts for all tests and benchmarks that we ran for the manuscript, in particular for assessing correctness and performance of grenedalf.

## Running the tests on grenedalf

The Python dependencies for running the tests are listed in `conda.yaml` and `pip.txt`, for the two package managers. Use these to install the packages needed to run the tests. To run the tests locally, make sure all Python dependencies are met. Then, the `execute_tests.py` script runs all 960 test cases, and the `evaluate.py` script checks that both implementations (grenedalf and the independent minimalistic Python implementation) yield the same results, as well as plots these against each other, for all estimators that we are interested in.

## Running the tests on PoPoolation

By default, we do not test PoPoolation here, as that's not the scope of this test. We however have implemented this for completeness as well; if you want to run the test for PoPoolation, find the line `with_popoolation = False` at the beginning of the `execute_test.py` script, and set it to `True`. Furthermore, download the source code of [PoPoolation](https://sourceforge.net/projects/popoolation/) and [PoPoolation2](https://sourceforge.net/projects/popoolation2/), and place both as sub-directories in `test/popoolation`, named `popoolation` and `popoolation2`, respectively.
