#!/usr/bin/env python3

# ================================================================================
#
# Imports
#
# ================================================================================
import argparse
import numpy.testing as npt
import pandas as pd


# ================================================================================
#
# Main
#
# ================================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare fextrema output")
    parser.add_argument(
        "-t", "--tol", help="Absolute tolerance", type=float, default=1e-15
    )
    parser.add_argument("-f", "--fname", help="Input file", type=str, required=True)
    parser.add_argument("-g", "--gold", help="Gold file", type=str, required=True)
    args = parser.parse_args()

    # Read the files
    cols = ["variables", "minimum_value", "maximum_value"]
    inpt = pd.read_csv(
        args.fname, delim_whitespace=True, skiprows=3, header=None, names=cols
    ).sort_values(by=["variables"])
    gold = pd.read_csv(
        args.gold, delim_whitespace=True, skiprows=3, header=None, names=cols
    ).sort_values(by=["variables"])

    # Compare
    npt.assert_allclose(inpt.minimum_value, gold.minimum_value, atol=args.tol)
    npt.assert_allclose(inpt.maximum_value, gold.maximum_value, atol=args.tol)
