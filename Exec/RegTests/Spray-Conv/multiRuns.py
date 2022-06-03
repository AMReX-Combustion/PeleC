#!/usr/bin/env python3

# A template script to lunch several times PeleC at different resolution
# in order to evaluate the convergence order
# Used in the marex regression_testing framework

# Usage:
#   ./multiRuns.py --test_name DummyTest --input_file PeleCInputFile

# Input:
#   * test_name: a TESTNAME that will be the name of a directory where the files go
#   * PeleCInputFile: the PeleC input file

# "Internal" user input
#   * bsizes : box sizes to test

# Head's up : 
#   * The PeleC executable is searched for in the current directory.

import sys
import os
import shutil
import argparse
import numpy as np

USAGE = """    
    A template script to launch several times PeleC.
"""

def multiRun(args):
    print(" Scripted runs for spray CI ")
    # User data
    bsizes = [32, 64]
    iter_1 = 5
    iter_2 = 9

    # Get the PeleC exec
    run_dir = os.getcwd()
    for f in os.listdir(run_dir):
        if ( f.startswith("PeleC") and f.endswith(".ex")):
               executable = f

    # Check the test name: current folder name is default
    if ( args.test_name == "None" ):
        args.test_name = "testfiles"

    # Check for the input file: first input.* is default
    if ( args.input_file == "None" ):
        for f in os.listdir(run_dir):
            if ( f.startswith("input") ):
                args.input_file = f
                break
    # Delete previous test directory if it exists
    for f in os.listdir(run_dir):
        if (f == args.test_name):
            os.system("rm -r {}".format(args.test_name))
    os.system("mkdir -p {}".format(args.test_name))
    for box in bsizes:
        print(" Running {}x{} box sizes".format(box,box))
        runtime_params = get_runtime_params(box, args)
        first_params = "max_step = {} ".format(iter_2)
        first_params += "amr.check_int = {} ".format(iter_1)
        first_params += "amr.checkpoint_files_output = 1 "
        first_params += "amr.plot_file = {}/{}_1_plt ".format(args.test_name, box)
        # Run first test to iter_2 and save a checkpoint at iter_1
        os.system("mpiexec -n 1 ./{} {} {} {}".format(executable, args.input_file, runtime_params, first_params))
        second_params = "max_step = {} ".format(iter_2)
        second_params += "amr.restart = {}/{}_chk000{} ".format(args.test_name, box, iter_1)
        second_params += "amr.checkpoint_files_output = 0 "
        second_params += "amr.plot_file = {}/{}_2_plt ".format(args.test_name, box)
        # Restart a case from checkpoint at iter_1
        os.system("mpiexec -n 1 ./{} {} {} {}".format(executable, args.input_file, runtime_params, second_params))
        # Now delete the checkpoint files and starting plot files
        print("Deleting unused plot and checkpoint files")
        os.system("rm -r {}/*plt0000".format(args.test_name, box))
        os.system("rm -r {}/{}_chk*".format(args.test_name, box))

def get_runtime_params(box, args):
    runtime_params = "amr.check_file={}/{}_chk ".format(args.test_name, box)
    runtime_params += "amr.plot_int=100 "
    runtime_params += "amr.initial_grid_file={}/gridfile_{}_1 ".format(args.grid_files_loc, box)
    runtime_params += "amr.regrid_file={}/gridfile_{}_2 ".format(args.grid_files_loc, box)
    runtime_params += "amr.max_grid_size={} ".format(box)
    runtime_params += "amr.blocking_factor={} ".format(box)
    return runtime_params

def parse_args(arg_string=None):
    parser = argparse.ArgumentParser(description=USAGE)

    parser.add_argument("--test_name", type=str, default="testfiles",
                        help="name of the test. Default = testfiles")

    parser.add_argument("--input_file", type=str, default="None",metavar="input.2d",
                        help="input file name. Default = first inputs.* in current directory")

    parser.add_argument("--grid_files_loc", type=str, default="two_d_gridfiles",
                        help="location of gridfiles. Default = two_d_gridfiles")

    if not arg_string is None:
        args, unknown = parser.parse_known_args(arg_string)
    else:
        args, unknown = parser.parse_known_args()

    return args

if __name__ == "__main__":
    args = parse_args(arg_string=sys.argv[1:])
    multiRun(args)
