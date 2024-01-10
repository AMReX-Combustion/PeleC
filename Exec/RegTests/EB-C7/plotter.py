import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

plt.rc("text", usetex=True)
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

labels = {
    "density": r"$\rho$",
    "pressure": r"$p$",
    "u": r"$u$",
    "v": r"$v$",
    "Temp": r"$T$",
}

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract data")
    parser.add_argument(
        "-f",
        "--fdirs",
        help="folders to post-process",
        type=str,
        required=True,
        nargs="+",
    )
    args = parser.parse_args()

    fields = ["density", "pressure", "u", "Temp"]

    exact = pd.read_csv("exact_soln", delim_whitespace=True)
    for field in fields:
        plt.figure(field)
        p = plt.plot(exact.x, exact[field], lw=1, color=cmap[-1], label=f"Exact")
        p[0].set_dashes(dashseq[-1])

    for k, fdir in enumerate(args.fdirs):
        df = pd.read_csv(os.path.join(fdir, "profiles.csv"))
        max_lev = df.max_level.iloc[0]
        dx = df.dxmax.iloc[0]
        for field in fields:
            plt.figure(field)
            p = plt.plot(
                df.xp,
                df[field],
                lw=2,
                color=cmap[k],
                label=f"$l = {max_lev}$, $\Delta x_0 = {dx}$",
            )
            p[0].set_dashes(dashseq[k])

    fname = "profiles.pdf"
    with PdfPages(fname) as pdf:
        for field in fields:
            plt.figure(field)
            ax = plt.gca()
            plt.xlabel(r"$x$", fontsize=22)
            plt.ylabel(labels[field], fontsize=22)
            plt.setp(ax.get_xmajorticklabels(), fontsize=16)
            plt.setp(ax.get_ymajorticklabels(), fontsize=16)
            plt.xlim(-0.5, 0.5)
            legend = ax.legend(loc="best")
            plt.tight_layout()
            pdf.savefig(dpi=300)
