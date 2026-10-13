# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# ParallelTestRunner runs each test_*.jl file in parallel, in a worker process.
# Each file is a module with a function for each test, and can also be run on
# its own. The methods directory has a file for each formulation.
#
# Positional arguments select the tests whose names start with them, for
# example
#
#     Pkg.test("PiecewiseLinearOpt"; test_args = ["methods/bivariate"])
#
# and --help lists the other options.
#
# The standard run tests every method, but combines the bivariate grid methods
# only with a few representative axis methods. The full run uses every
# combination; start it with
#
#     Pkg.test("PiecewiseLinearOpt"; test_args = ["--full"])
#
# or by setting the environment variable PLO_TEST_FULL=true.

import PiecewiseLinearOpt

using ParallelTestRunner

args = parse_args(ARGS; custom = ["full"])
if args.custom["full"] !== nothing
    # The workers inherit the environment
    ENV["PLO_TEST_FULL"] = "true"
end

testsuite = find_tests(@__DIR__)
filter!(((name, _),) -> startswith(basename(name), "test_"), testsuite)
filter_tests!(testsuite, args)

runtests(PiecewiseLinearOpt, args; testsuite)
