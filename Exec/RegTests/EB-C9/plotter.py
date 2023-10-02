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


def parse_ic(fname):
    """
    Parse the file written by PeleC to understand the initial condition

    Returns a dictionary for easy access
    """

    # Read into dataframe
    df = pd.read_csv(fname)
    df.rename(columns=lambda x: x.strip(), inplace=True)

    # convert to dictionary for easier access
    return df.to_dict("records")[0]


def theory_ooa(order, res, orig):
    return orig * (res[0] / res) ** order


def rho_exact(ics, x):
    return ics["rho"] + ics["alpha"] * np.exp(
        -(((x - ics["cs"] * tf) / ics["sigma"]) ** 2)
    )


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

    ics = parse_ic(os.path.join(args.fdirs[0], "ic.txt"))

    errors = np.zeros((2, len(args.fdirs)))
    for k, fdir in enumerate(args.fdirs):
        df = pd.read_csv(os.path.join(fdir, "profiles.csv"))
        df["entropy"] = df.pressure / (df.density ** ics["gamma"])
        final = df[df.time == df.time.max()].reset_index()
        tf = df.time.max()

        res = int(fdir)
        errors[0, k] = res
        errors[1, k] = np.sqrt(
            np.sum((final.density - rho_exact(ics, final.x)) ** 2) / res
        )

        plt.figure("rho")
        if k == 0:
            x = np.linspace(-50, 50, 1000)
            rho_e = rho_exact(ics, x)
            plt.plot(x, rho_e, lw=2, color=cmap[-1], label=r"Exact")

        p = plt.plot(
            final.x, final.density, lw=2, color=cmap[k], label=f"$n_x = {res}$"
        )
        p[0].set_dashes(dashseq[k])

        plt.figure("entropy")
        p = plt.plot(
            final.x,
            final.entropy / (ics["p"] / (ics["rho"] ** ics["gamma"])),
            lw=2,
            color=cmap[k],
            label=f"$n_x = {res}$",
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

        plt.figure("entropy")
        ax = plt.gca()
        plt.xlabel(r"$x~[cm]$", fontsize=22)
        plt.ylabel(r"$(p / \rho^\gamma) / (p_0 / \rho_0^\gamma)$", fontsize=22)
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
