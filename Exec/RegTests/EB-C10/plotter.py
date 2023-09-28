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


def eval_u_exact(ics, r, factor):
    return factor * ics["G"] / (4.0 * ics["mu"]) * (ics["radius"] ** 2 - r**2)


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

    ics = parse_ic(os.path.join(args.fdirs[0], "ic.txt"))
    r = np.linspace(-ics["radius"], ics["radius"], 1000)
    factor = 0.993
    u_exact = eval_u_exact(ics, r, factor)

    plt.figure("u")
    plt.plot(
        r / ics["radius"], u_exact / ics["umax"], lw=2, color=cmap[-1], label=r"Exact"
    )

    errors = np.zeros((2, len(args.fdirs)))
    for k, fdir in enumerate(args.fdirs):
        df = pd.read_csv(os.path.join(fdir, "profiles.csv"))
        df = df[df.xloc == 6.0].reset_index()

        res = int(fdir)
        errors[0, k] = 2 * res
        errors[1, k] = np.sqrt(
            np.sum(
                (
                    df.x_velocity / ics["umax"]
                    - eval_u_exact(ics, df.y, factor) / ics["umax"]
                )
                ** 2
            )
            / (2 * res)
        )

        plt.figure("u")
        p = plt.plot(
            df.y / ics["radius"],
            df.x_velocity / ics["umax"],
            lw=2,
            color=cmap[k],
            label=f"$n_r = {res}$",
        )
        p[0].set_dashes(dashseq[k])

        df = pd.read_csv(os.path.join(fdir, "history.csv"))
        plt.figure("ke")
        plt.plot(
            df.time / (ics["L"] / ics["umax"]),
            df.rho_ke
            / (
                0.5
                * ics["rho"]
                * (0.5 * ics["umax"]) ** 2
                * ics["L"]
                * np.pi
                * ics["radius"] ** 2
            ),
            lw=2,
            color=cmap[k],
            label=f"$n_r = {res}$",
        )

    plt.figure("error")
    plt.loglog(
        errors[0, :],
        errors[1, :],
        lw=2,
        color=cmap[0],
        marker=markertype[0],
        label=r"Sim.",
    )

    p1 = theory_ooa(1, errors[0, :], errors[1, 0])
    p2 = theory_ooa(2, errors[0, :], errors[1, 0])
    plt.loglog(errors[0, :], p1, lw=2, color=cmap[-1], label=f"$p=1$")
    plt.loglog(errors[0, :], p2, lw=2, color=cmap[-2], label=f"$p=2$")

    print("Estimated order of the error:")
    print(np.log(errors[1, :-1] / errors[1, 1:]) / np.log(2))

    fname = "profiles.pdf"
    with PdfPages(fname) as pdf:
        plt.figure("u")
        ax = plt.gca()
        plt.xlabel(r"$r / R$", fontsize=22)
        plt.ylabel(r"$u / u_m$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("error")
        ax = plt.gca()
        plt.xlabel(r"$2 n_r$", fontsize=22)
        plt.ylabel(r"$\epsilon$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("ke")
        ax = plt.gca()
        plt.xlabel(r"$t / \tau$", fontsize=22)
        plt.ylabel(r"$K / K_e$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)
