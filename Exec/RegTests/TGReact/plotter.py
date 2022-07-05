#!/usr/bin/env python3

import argparse
import os
import re
import glob as glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd


plt.rc("text", usetex=True)
plt.rc("figure", max_open_warning=100)
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


def parse_ic(fname):
    """
    Parse the file written by PeleC to understand the initial condition

    Returns a dictionary for easy acces
    """

    # Read into dataframe
    df = pd.read_csv(fname)
    df.rename(columns=lambda x: x.strip(), inplace=True)

    # convert to dictionary for easier access
    return df.to_dict("records")[0]


def plot_2d_refdata(ics):

    time = 10 * ics["tau"]
    decay = np.exp(
        -8 * (np.pi**2) * (ics["mu"] / ics["rho0"]) * time / (ics["L"] ** 2)
    )

    x = np.linspace(0, ics["L"], 100)
    y = 0.5 * ics["L"]
    ux = (
        ics["u0"]
        * np.sin(ics["omega_x"] * x / ics["L"])
        * np.cos(ics["omega_y"] * y / ics["L"])
        * decay
    )
    plt.figure("ux")
    plt.plot(x / ics["L"], ux / ics["u0"], ls="-", color=cmap[-1], lw=1, label="Ref.")

    x = 0.5 * ics["L"]
    y = np.linspace(0, ics["L"], 100)
    uy = (
        -ics["u0"]
        * np.cos(ics["omega_x"] * x / ics["L"])
        * np.sin(ics["omega_y"] * y / ics["L"])
        * decay
    )
    plt.figure("uy")
    plt.plot(y / ics["L"], uy / ics["u0"], ls="-", color=cmap[-1], lw=1, label="Ref.")


def plot_cold3d_refdata(refdir, ics):
    m2cm = 100

    KE_ref = pd.read_csv(os.path.join(refdir, "KE.csv"))
    plt.figure("KE")
    plt.plot(KE_ref.t, KE_ref.KE, ls="-", color=cmap[-1], lw=1, label="Ref.")

    epsilon_ref = pd.read_csv(os.path.join(refdir, "epsilon.csv"))
    plt.figure("epsilon")
    plt.plot(
        epsilon_ref.t, epsilon_ref.epsilon, ls="-", color=cmap[-1], lw=1, label="Ref."
    )

    ux_ref = pd.read_csv(os.path.join(refdir, "ux.csv"))
    ux_ref.x *= m2cm
    ux_ref.ux *= m2cm
    plt.figure("ux")
    plt.plot(
        ux_ref.x / ics["L"],
        ux_ref.ux / ics["u0"],
        ls="-",
        color=cmap[-1],
        lw=1,
        label="Ref.",
    )

    uy_ref = pd.read_csv(os.path.join(refdir, "uy.csv"))
    uy_ref.y *= m2cm
    uy_ref.uy *= m2cm
    plt.figure("uy")
    plt.plot(
        uy_ref.y / ics["L"],
        uy_ref.uy / ics["u0"],
        ls="-",
        color=cmap[-1],
        lw=1,
        label="Ref.",
    )

    vortz_ref = pd.read_csv(os.path.join(refdir, "vortz.csv"))
    vortz_ref.x *= m2cm
    plt.figure("vortz")
    plt.plot(
        vortz_ref.x / ics["L"],
        vortz_ref.vortz,
        ls="-",
        color=cmap[-1],
        lw=1,
        label="Ref.",
    )


def plot_nonreacting3d_refdata(refdir, ics):
    ms2s = 1e-3
    u0 = 4  # m/s

    ux = pd.read_csv(os.path.join(refdir, "ux_nonreacting.csv"))
    plt.figure("ux")
    plt.plot(ux.x / ics["L"], ux.ux / u0, ls="-", color=cmap[-1], lw=1, label="Ref.")

    uy = pd.read_csv(os.path.join(refdir, "uy_nonreacting.csv"))
    plt.figure("uy")
    plt.plot(uy.y / ics["L"], uy.uy / u0, ls="-", color=cmap[-1], lw=1, label="Ref.")

    temp = pd.read_csv(os.path.join(refdir, "temp_nonreacting.csv"))
    plt.figure("temp")
    plt.plot(temp.y / ics["L"], temp.temp, ls="-", color=cmap[-1], lw=1, label="Ref.")

    YH2 = pd.read_csv(os.path.join(refdir, "YH2_nonreacting.csv"))
    plt.figure("YH2")
    plt.plot(YH2.y / ics["L"], YH2.YH2, ls="-", color=cmap[-1], lw=1, label="Ref.")

    YO2 = pd.read_csv(os.path.join(refdir, "YO2_nonreacting.csv"))
    plt.figure("YO2")
    plt.plot(YO2.y / ics["L"], YO2.YO2, ls="-", color=cmap[-1], lw=1, label="Ref.")

    maxtemp = pd.read_csv(os.path.join(refdir, "maxtemp_nonreacting.csv"))
    maxtemp.t *= ms2s
    plt.figure("maxtemp")
    plt.plot(
        maxtemp.t / ics["tau"], maxtemp.temp, ls="-", color=cmap[-1], lw=1, label="Ref."
    )

    temp = pd.read_csv(os.path.join(refdir, "temp_init_nonreacting.csv"))
    plt.figure("temp_init")
    plt.plot(temp.x / ics["L"], temp.temp, ls="-", color=cmap[-1], lw=1, label="Ref.")

    YH2 = pd.read_csv(os.path.join(refdir, "YH2_init_nonreacting.csv"))
    plt.figure("YH2_init")
    plt.plot(YH2.x / ics["L"], YH2.YH2, ls="-", color=cmap[-1], lw=1, label="Ref.")

    YO2 = pd.read_csv(os.path.join(refdir, "YO2_init_nonreacting.csv"))
    plt.figure("YO2_init")
    plt.plot(YO2.x / ics["L"], YO2.YO2, ls="-", color=cmap[-1], lw=1, label="Ref.")


def plot_reacting3d_refdata(refdir, ics):
    ms2s = 1e-3
    u0 = 4  # m/s

    ux = pd.read_csv(os.path.join(refdir, "ux_reacting.csv"))
    plt.figure("ux")
    plt.plot(ux.x / ics["L"], ux.ux / u0, ls="-", color=cmap[-1], lw=1, label="Ref.")

    uy = pd.read_csv(os.path.join(refdir, "uy_reacting.csv"))
    plt.figure("uy")
    plt.plot(uy.y / ics["L"], uy.uy / u0, ls="-", color=cmap[-1], lw=1, label="Ref.")

    temp = pd.read_csv(os.path.join(refdir, "temp_reacting.csv"))
    plt.figure("temp")
    plt.plot(temp.y / ics["L"], temp.temp, ls="-", color=cmap[-1], lw=1, label="Ref.")

    YH2 = pd.read_csv(os.path.join(refdir, "YH2_reacting.csv"))
    plt.figure("YH2")
    plt.plot(YH2.y / ics["L"], YH2.YH2, ls="-", color=cmap[-1], lw=1, label="Ref.")

    YO2 = pd.read_csv(os.path.join(refdir, "YO2_reacting.csv"))
    plt.figure("YO2")
    plt.plot(YO2.y / ics["L"], YO2.YO2, ls="-", color=cmap[-1], lw=1, label="Ref.")

    YOH = pd.read_csv(os.path.join(refdir, "YOH_reacting.csv"))
    plt.figure("YOH")
    plt.plot(YOH.y / ics["L"], YOH.YOH, ls="-", color=cmap[-1], lw=1, label="Ref.")

    HR = pd.read_csv(os.path.join(refdir, "HR_reacting.csv"))
    plt.figure("HR")
    plt.plot(HR.y / ics["L"], HR.HR, ls="-", color=cmap[-1], lw=1, label="Ref.")

    maxtemp = pd.read_csv(os.path.join(refdir, "maxtemp_reacting.csv"))
    maxtemp.t *= ms2s
    plt.figure("maxtemp")
    plt.plot(
        maxtemp.t / ics["tau"], maxtemp.temp, ls="-", color=cmap[-1], lw=1, label="Ref."
    )

    temp = pd.read_csv(os.path.join(refdir, "temp_init_reacting.csv"))
    plt.figure("temp_init")
    plt.plot(temp.x / ics["L"], temp.temp, ls="-", color=cmap[-1], lw=1, label="Ref.")

    YH2 = pd.read_csv(os.path.join(refdir, "YH2_init_reacting.csv"))
    plt.figure("YH2_init")
    plt.plot(YH2.x / ics["L"], YH2.YH2, ls="-", color=cmap[-1], lw=1, label="Ref.")

    YO2 = pd.read_csv(os.path.join(refdir, "YO2_init_reacting.csv"))
    plt.figure("YO2_init")
    plt.plot(YO2.x / ics["L"], YO2.YO2, ls="-", color=cmap[-1], lw=1, label="Ref.")


if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(description="A simple plot tool")
    parser.add_argument(
        "-f",
        "--folders",
        nargs="+",
        help="Folders where files are stored",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-c",
        "--case",
        help="Type of problem",
        type=str,
        default="cold3d",
        choices=["2d", "cold3d", "nonreacting3d", "reacting3d"],
    )
    args = parser.parse_args()

    # Setup
    icname = "ic.txt"
    datname = "datlog"
    extraname = "extralog"
    plot_time = 0.0  # (in PeleC units)
    visit_time_fields = ["maxtemp", "Sij2", "magvel2"]
    profile_fields = {"z_vorticity": ["vortz", 1.0], "Temp": ["temp", 1.0]}

    # Loop on folders
    for i, folder in enumerate(args.folders):

        # Setup
        fdir = os.path.abspath(folder)
        datadir = os.path.join(fdir, "data")
        ics = parse_ic(os.path.join(fdir, icname))
        resolution = int(folder)
        dx = ics["L"] / resolution
        volume = ics["L"] ** 3
        profile_fields["x_velocity"] = ["ux", ics["u0"]]
        profile_fields["y_velocity"] = ["uy", ics["u0"]]

        # Set case specific stuff
        if args.case == "2d":
            plot_time = 2.0
        elif args.case == "cold3d":
            plot_time = 12.11
            ics["tau"] = 1.0
        elif args.case in ["nonreacting3d", "reacting3d"]:
            plot_time = 0.5 * 1e-3
            ics["tau"] = 1.0
            profile_fields["Y_lp_H2_rp_"] = ["YH2", 1.0]
            profile_fields["Y_lp_O2_rp_"] = ["YO2", 1.0]
            profile_fields["init_Temp"] = ["temp_init", 1.0]
            profile_fields["init_Y_lp_H2_rp_"] = ["YH2_init", 1.0]
            profile_fields["init_Y_lp_O2_rp_"] = ["YO2_init", 1.0]
            if args.case == "reacting3d":
                profile_fields["Y_lp_OH_rp_"] = ["YOH", 1.0]
                profile_fields["heat_release"] = ["HR", -10.0]

        # Read time history files
        datdf = pd.read_csv(os.path.join(fdir, datname), delim_whitespace=True)
        extradf = pd.read_csv(os.path.join(fdir, extraname), delim_whitespace=True)
        tdf = pd.concat(
            [datdf.set_index("time"), extradf.set_index("time")], axis=1, join="inner"
        ).reset_index()
        tdf["KE"] = tdf.rho_K / tdf.rho_K.iloc[0]
        tdf["t"] = tdf.time / ics["tau"]
        tdf.sort_values(by=["t"], inplace=True)

        # Read visit time history files
        vtdf = pd.concat(
            [
                pd.read_csv(
                    os.path.join(datadir, f"{field}.curve"),
                    delim_whitespace=True,
                    header=None,
                    names=["t", field],
                    comment="#",
                ).set_index("t")
                for field in visit_time_fields
            ],
            axis=1,
            join="inner",
        ).reset_index()
        vtdf["KE"] = 0.5 * vtdf.magvel2 / volume
        vtdf["KE"] /= vtdf.KE.iloc[0]
        vtdf["t"] /= ics["tau"]
        vtdf.sort_values(by=["t"], inplace=True)

        # Read visit profile files
        profiles = {}
        for field in profile_fields.keys():
            fnames = glob.glob(os.path.join(datadir, f"{field}*.curve"))
            lst = []
            for j, fname in enumerate(fnames):
                tmp = pd.read_csv(
                    fname,
                    delim_whitespace=True,
                    header=None,
                    names=["pos", profile_fields[field][0]],
                    comment="#",
                )
                # nondimensionalize and offset to cell center
                tmp["pos"] = (tmp["pos"] + 0.5 * dx) / ics["L"]
                tmp[profile_fields[field][0]] /= profile_fields[field][1]
                tmp["step"] = int(os.path.splitext(fname)[0].split("_")[-2])
                tmp["time"] = float(os.path.splitext(fname)[0].split("_")[-1])
                lst.append(tmp)

            profiles[profile_fields[field][0]] = pd.concat(lst).sort_values(
                by=["step", "pos"]
            )

        # Plot for time quantities
        for field in ["KE", "maxtemp"]:
            plt.figure(field)
            p = plt.plot(tdf.t, tdf[field], ls="-", color=cmap[i], lw=2)
            p[0].set_dashes(dashseq[i])

        # Plot for visit time quantities (just epsilon actually)
        if args.case in ["2d", "cold3d"]:
            field = "epsilon"
            vtdf[field] = 2.0 * ics["mu"] / ics["rho0"] * vtdf.Sij2
            vtdf[field] /= vtdf.epsilon.max()

            plt.figure(field)
            p = plt.plot(vtdf.t, vtdf[field], ls="-", color=cmap[i], lw=2)
            p[0].set_dashes(dashseq[i])

        # Plot all the profiles at a specific time
        plot_step = profiles["ux"].step.unique()[
            np.argmin(np.fabs(profiles["ux"].time.unique() - plot_time))
        ]
        for key, df in profiles.items():
            plt.figure(key)
            sdf = df[df.step == (plot_step if "init" not in key else 0)]
            p = plt.plot(sdf.pos, sdf[key], ls="-", color=cmap[i], lw=2)
            p[0].set_dashes(dashseq[i])

    # Plot reference data
    refdir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "refdata")
    if args.case == "2d":
        plot_2d_refdata(ics)
    elif args.case == "cold3d":
        plot_cold3d_refdata(refdir, ics)
    elif args.case == "nonreacting3d":
        plot_nonreacting3d_refdata(refdir, ics)
    elif args.case == "reacting3d":
        plot_reacting3d_refdata(refdir, ics)

    # Save plots
    fname = "plots.pdf"
    with PdfPages(fname) as pdf:

        plt.figure("KE")
        ax = plt.gca()
        plt.xlabel(r"$t / \tau$", fontsize=22)
        plt.ylabel(r"$K / K_0$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("epsilon")
        ax = plt.gca()
        plt.xlabel(r"$t / \tau$", fontsize=22)
        plt.ylabel(r"$\epsilon / \epsilon_{\max}$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("maxtemp")
        ax = plt.gca()
        plt.xlabel(r"$t / \tau$", fontsize=22)
        plt.ylabel(r"$T_{\max}~[K]$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("ux")
        ax = plt.gca()
        plt.xlabel(r"$x / L$", fontsize=22)
        plt.ylabel(r"$u_x / u_0$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("uy")
        ax = plt.gca()
        plt.xlabel(r"$y / L$", fontsize=22)
        plt.ylabel(r"$u_y / u_0$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("vortz")
        ax = plt.gca()
        plt.xlabel(r"$x / L$", fontsize=22)
        plt.ylabel(r"$\omega_z~[1/s]$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("temp")
        ax = plt.gca()
        plt.xlabel(r"$y / L$", fontsize=22)
        plt.ylabel(r"$T~[K]$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("YH2")
        ax = plt.gca()
        plt.xlabel(r"$y / L$", fontsize=22)
        plt.ylabel(r"$Y_{H_2}$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("YO2")
        ax = plt.gca()
        plt.xlabel(r"$y / L$", fontsize=22)
        plt.ylabel(r"$Y_{O_2}$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("YOH")
        ax = plt.gca()
        plt.xlabel(r"$y / L$", fontsize=22)
        plt.ylabel(r"$Y_{OH}$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("HR")
        ax = plt.gca()
        plt.xlabel(r"$y / L$", fontsize=22)
        plt.ylabel(r"$HR~[W/m^3]$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("temp_init")
        ax = plt.gca()
        plt.xlabel(r"$x / L$", fontsize=22)
        plt.ylabel(r"$T~[K]$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("YH2_init")
        ax = plt.gca()
        plt.xlabel(r"$x / L$", fontsize=22)
        plt.ylabel(r"$Y_{H_2}$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("YO2_init")
        ax = plt.gca()
        plt.xlabel(r"$x / L$", fontsize=22)
        plt.ylabel(r"$Y_{O_2}$", fontsize=22)
        plt.setp(ax.get_xmajorticklabels(), fontsize=16)
        plt.setp(ax.get_ymajorticklabels(), fontsize=16)
        plt.tight_layout()
        pdf.savefig(dpi=300)
