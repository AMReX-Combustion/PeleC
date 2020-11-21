#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ========================================================================
#
# Imports
#
# ========================================================================
import os
import re
import numpy as np
import numpy.testing as npt
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import unittest

# ========================================================================
#
# Some defaults variables
#
# ========================================================================
# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', serif='Times')
cmap_med = [
    "#F15A60",
    "#7AC36A",
    "#5A9BD4",
    "#FAA75B",
    "#9E67AB",
    "#CE7058",
    "#D77FB4",
    "#737373",
]
cmap = [
    "#EE2E2F",
    "#008C48",
    "#185AA9",
    "#F47D23",
    "#662C91",
    "#A21D21",
    "#B43894",
    "#010202",
]
dashseq = [
    (None, None),
    [10, 5],
    [10, 4, 3, 4],
    [3, 3],
    [10, 4, 3, 4, 3, 4],
    [3, 3],
    [3, 3],
]
markertype = ["s", "d", "o", "p", "h"]


# ========================================================================
#
# Function definitions
#
# ========================================================================
def load_pelec_error(fdir, theory_order):
    """Load the error for each resolution"""
    lst = []
    resolutions = sorted(
        [
            int(f)
            for f in os.listdir(fdir)
            if os.path.isdir(os.path.join(fdir, f)) and re.match("^[0-9]+$", f)
        ],
        key=int,
    )
    resdirs = [os.path.join(fdir, str(res)) for res in resolutions]

    for k, (res, resdir) in enumerate(zip(resolutions, resdirs)):

        fname = os.path.join(resdir, "mmslog")
        df = pd.read_csv(fname, delim_whitespace=True)

        idx = -1
        print(
            "Loading {0:d} at t = {1:e} (step = {2:d})".format(
                res, df["time"].iloc[idx], df.index[idx]
            )
        )
        lst.append(
            [
                res,
                1.0 / res,
                df["rho_mms_err"].iloc[idx],
                df["u_mms_err"].iloc[idx],
                df["v_mms_err"].iloc[idx],
                df["w_mms_err"].iloc[idx],
                df["p_mms_err"].iloc[idx],
            ]
        )

    edf = pd.DataFrame(
        lst,
        columns=[
            "resolution",
            "dx",
            "rho_mms_err",
            "u_mms_err",
            "v_mms_err",
            "w_mms_err",
            "p_mms_err",
        ],
    )

    # Theoretical error
    idx = 1
    edf["rho_theory"] = (
        edf["rho_mms_err"].iloc[idx]
        * (edf["resolution"].iloc[idx] / edf["resolution"]) ** theory_order
    )
    edf["u_theory"] = (
        edf["u_mms_err"].iloc[idx]
        * (edf["resolution"].iloc[idx] / edf["resolution"]) ** theory_order
    )
    edf["v_theory"] = (
        edf["v_mms_err"].iloc[idx]
        * (edf["resolution"].iloc[idx] / edf["resolution"]) ** theory_order
    )
    edf["w_theory"] = (
        edf["w_mms_err"].iloc[idx]
        * (edf["resolution"].iloc[idx] / edf["resolution"]) ** theory_order
    )
    edf["p_theory"] = (
        edf["p_mms_err"].iloc[idx]
        * (edf["resolution"].iloc[idx] / edf["resolution"]) ** theory_order
    )

    return edf.loc[:, (edf != 0).any(axis=0)]


def calculate_ooa(edf):
    """Calculate the order of accuracy given an error dataframe."""

    sfx_mms = "_mms_err"
    fields = [re.sub(sfx_mms, "", col) for col in edf.columns if col.endswith(sfx_mms)]

    columns = []
    data = np.zeros((len(edf["resolution"]) - 1, len(fields)))
    for k, field in enumerate(fields):
        columns.append(field + "_ooa")
        data[:, k] = -np.diff(np.log(edf[field + sfx_mms])) / np.diff(
            np.log(edf["resolution"])
        )
    ooa = pd.DataFrame(data, columns=columns)

    return ooa.dropna(axis=1, how="any")


def plot_errors(fdir, edf):
    """Plot the error dataframe."""

    sfx_mms = "_mms_err"
    fields = [re.sub(sfx_mms, "", col) for col in edf.columns if col.endswith(sfx_mms)]

    plt.close("all")
    for k, field in enumerate(fields):

        plt.figure(k)
        p = plt.loglog(
            edf["resolution"],
            edf[field + sfx_mms],
            ls="-",
            lw=2,
            color=cmap[0],
            marker=markertype[0],
            mec=cmap[0],
            mfc=cmap[0],
            ms=10,
            label="Pele",
        )
        p = plt.loglog(
            edf["resolution"],
            edf[field + "_theory"],
            ls="-",
            lw=2,
            color=cmap[-1],
            label="2nd order",
        )

        # Format the plots
        ax = plt.gca()
        plt.xlabel(r"$N$", fontsize=22, fontweight="bold")
        plt.ylabel(
            r"$e_{0:s}$".format("{" + field + "}"), fontsize=22, fontweight="bold"
        )
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        legend = ax.legend(loc="best")
        plt.tight_layout()
        plt.savefig(os.path.join(fdir, field + "_error.png"), format="png")


# ========================================================================
#
# Test definitions
#
# ========================================================================
class MMSTestCase(unittest.TestCase):
    """Tests for the order of accuracy in Pele."""

    def setUp(self):

        self.theory_order = 2.0
        self.tol = 2.5e-1

    def test_second_order(self):
        """Is this test second order accurate?"""

        # Load the data
        fdir = os.path.abspath(".")
        edf = load_pelec_error(fdir, self.theory_order)
        ooa = calculate_ooa(edf)

        # Plot the errors
        plot_errors(fdir, edf)

        # Test against theoretical OOA
        npt.assert_allclose(np.array(ooa.iloc[-1]) > self.theory_order - self.tol, True)


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == "__main__":
    unittest.main()
