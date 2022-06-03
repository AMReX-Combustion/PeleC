#!/usr/bin/env python3

# Template post-processing script for PeleC/spray regression testing
# Must be used after multirun.py script
# Input are limited by the regression framework.

# Usage:
#   ./compareOutput.py --test_name DummyTest

# Input:
#   * --test_name: a TESTNAME that will looked for during the postprocessing
#   * --max_error: the maximum allow able error

# "Internal" user input
#   * vars : a list of the variables of interest

# Output:
#  If the error in vars is higher than the max_error, an error statement is returned

# Head's up :
#   - The script will get a copy of the post-processing program (if not already there) in the testing folder. The name of this folder is assumed to be the TESTNAME.
#   - The plt files naming convention is: ${TESTNAME}/64_1_plt***** where 64 is the box size and 1 is the first test
#   - Errors are parsed from the screen output of the standard fcompare. Beware of any change of these programs.

import sys
import os
import fnmatch
import shutil
import argparse
import numpy as np

USAGE = """
    Template post-processing script for PeleC error analysis
"""

def compare(args):

    # User data
    vars=["density", "rho_NC10H22", "xmom", "Temp", "rho_E", "spray_vol"]

    # Get a local copy of post-processing executable
    run_dir = os.getcwd()
    fcomp_exe = "None"
    for f in os.listdir(run_dir):
        if ( f.startswith("fcompare") and f.endswith(".ex")):
               fcomp_exe = f
    if (fcomp_exe == "None"):
        errorStatement="fcompare executable not found"
        raise ValueError(errorStatement)
    # Check the test name: current folder name is default
    if ( args.test_name == "None" ):
        args.test_name = "testfiles"

    # Run the postprocessing
    pltfiles = ['', '', '', '']
    # Put all plot files in list
    # Note: this assumes only the final plot files are there and plt00000 has been removed
    for f in os.listdir(args.test_name):
        if ( not fnmatch.fnmatch(f, '*old*')):
            if (f.startswith("32_1_plt")):
                pltfiles[0] = f
            elif (f.startswith("32_2_plt")):
                pltfiles[1] = f
            elif (f.startswith("64_1_plt")):
                pltfiles[2] = f
            elif (f.startswith("64_2_plt")):
                pltfiles[3] = f
    # We have 3 comparisons to make,
    # 32_1 to 64_1, 32_1 to 64_2, and 64_1 to 64_2
    numcomps = 3
    comp1 = [0, 0, 2]
    comp2 = [2, 3, 3]
    maxerror = args.max_error
    errfail = 0
    for comp in range(numcomps):
        indx1 = comp1[comp]
        indx2 = comp2[comp]
        plt1 = args.test_name + "/" + pltfiles[indx1]
        plt2 = args.test_name + "/" + pltfiles[indx2]
        print("Comparing the L2 norm of the error between " + plt1 + " and " + plt2)
        outfile = "error_{}.analysis.out".format(comp)
        os.system("./{} -n 2 -a {} {} > {}".format(os.path.basename(fcomp_exe), plt1, plt2, outfile))
        # Extract errors on each variable
        curlev = 0
        varcount = 0
        with open(outfile) as fp:
            for i, line in enumerate(fp):
                if (i >= 2):
                    var = line.split()[0]
                    if (var == "level"):
                        curlev = int(line.split()[2])
                    for v in range(len(vars)):
                        if ( var == vars[v] ):
                            if (curlev == 0):
                                varcount += 1
                            cer = float(line.split()[2])
                            if (cer > maxerror):
                                errfail += 1
                                perr = "Error in {} on level {}: {} > {}".format(var,curlev,cer,maxerror)
                                print(perr)
        if (varcount != len(vars)):
            error = "Not all variables were found in plot file"
            raise ValueError(error)
    if (errfail > 0):
        errorStatement = "Failure detected in spray test"
        raise ValueError(errorStatement)
    else:
        print("Spray tests passed")

def parse_args(arg_string=None):
    parser = argparse.ArgumentParser(description=USAGE)

    parser.add_argument("--test_name", type=str, default="testfiles",
                        help="name of the test. Default = testfiles")

    parser.add_argument("--max_error", type=float, default=1.E-12,
                        help="max error for tests.")

    if not arg_string is None:
        args, unknown = parser.parse_known_args(arg_string)
    else:
        args, unknown = parser.parse_known_args()

    return args

if __name__ == "__main__":
    args = parse_args(arg_string=sys.argv[1:])
    compare(args)
