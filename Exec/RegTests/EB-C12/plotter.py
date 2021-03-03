import os
import glob
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


def theory_ooa(order, res, orig):
    return orig * (res[0] / res) ** order


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

    x = np.linspace(-1, 1, 1000)
    rho_exact = 1 + 0.2 * np.sin(np.pi * x)
    plt.figure("rho")
    plt.plot(x, rho_exact, lw=2, color=cmap[-1], label=r"Exact")

    errors = np.zeros((2, len(args.fdirs)))
    for k, fdir in enumerate(args.fdirs):
        df = pd.read_csv(os.path.join(fdir, "profiles.csv"))
        init = df[df.time == df.time.min()].reset_index()
        final = df[df.time == df.time.max()].reset_index()

        res = int(fdir) * 10
        errors[0, k] = res
        errors[1, k] = np.sqrt(np.sum((final.density - init.density) ** 2) / res)

        plt.figure("rho")
        p = plt.plot(
            final.x, final.density, lw=2, color=cmap[k], label=f"$n_x = {res}$"
        )
        p[0].set_dashes(dashseq[k])

    plt.figure("error")
    plt.loglog(
        errors[0, :],
        errors[1, :],
        lw=2,
        color=cmap[0],
        marker=markertype[0],
        label=r"Sim.",
    )

    p2 = theory_ooa(2, errors[0, :], errors[1, 0])
    plt.loglog(errors[0, :], p2, lw=2, color=cmap[-1], label=f"$p=2$")

    print("Estimated order of the error:")
    print(np.log(errors[1, :-1] / errors[1, 1:]) / np.log(2))

    fname = "profiles.pdf"
    with PdfPages(fname) as pdf:
        plt.figure("rho")
        ax = plt.gca()
        plt.xlabel(r"$x~[cm]$", fontsize=22)
        plt.ylabel(r"$\rho$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("error")
        ax = plt.gca()
        plt.xlabel(r"$n_x$", fontsize=22)
        plt.ylabel(r"$\epsilon$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)
