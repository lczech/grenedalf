# Test

This directory contains test cases that are mostly meant for internal testing. However, feel free to have a look if you are interested. The created test data can for instance be useful as a standardized test for evaluating other tools.

We here run grenedalf as well as our minimalistic independent implementation of the equations in Python, and compare the two to each other. As they both yield the same results, we gain some confidence that the equations are correctly implemented. The tests span a range of pool sizes, read depths, allele frequencies, and window sizes. However, this is a minimalistic test with clean data, i.e., we assume sufficient read depth and no missing data in the test. Hence, if you encounter any issues on more realistic data, please report the [issue](https://github.com/lczech/grenedalf/issues).

The Python dependencies for running the tests are listed in `conda.yaml` and `pip.txt`, for the two package manages. Use these to install the packages needed to run the tests. We also run these tests in our [CI Tests](https://github.com/lczech/grenedalf/actions) in GitHub Actions.

To run the tests locally, make sure all Python dependencies are met. Then, the `execute_tests.py` script runs all 960 test cases, and the `evaluate.py` script checks that both implementations (grenedalf and the independent minimalistic Python implementation) yield the same results, as well as plots these against each other, for all estimators that we are interested in.

By default, we do not test PoPoolation here, as that's not the scope of this test. We however have implemented this for completeness as well; if you want to run the test for PoPoolation, find the commented line containing `run_popoolation` at the very end of the `execute_test.py` script, and un-comment it. Furthermore, download the source code of [PoPoolation](https://sourceforge.net/projects/popoolation/) and [PoPoolation2](https://sourceforge.net/projects/popoolation2/), and place both as sub-directories in `test/popoolation`, named `popoolation` and `popoolation2`, respectively.

The test scripts are also published in the supporting grenedalf manuscript repository, see [here](https://github.com/lczech/grenedalf-paper/tree/master/eval-independent-test). That repository contains the scripts for all tests and benchmarks that we ran for the manuscript, in particular for assessing correctness and performance of grenedalf.
