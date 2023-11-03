#!/bin/bash

# Run tests
R -e "x <- shiny::runTests(assert = FALSE); writeLines(as.character(all(x[[2]])), 'test_result.txt')"

# Read test results from file
test_result=$(cat test_result.txt)

# return test result as an output
#echo ::set-output name=test_result::$test_result
echo "{test_result}={$test_result}" >> $env:GITHUB_OUTPUT # return test result as an output

