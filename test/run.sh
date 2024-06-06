#!/bin/bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

echo "Running execute_tests.py"
${SCRIPT_DIR}/execute_tests.py
if [ $? -ne 0 ]; then
    echo "FAIL"
    exit 1
fi

echo "Running evaluate.py"
${SCRIPT_DIR}/evaluate.py
if [ $? -ne 0 ]; then
    echo "FAIL"
    exit 1
fi
