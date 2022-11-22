#!/usr/bin/env python3

# ========================================================================
#
# Imports
#
# ========================================================================
import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd


# ========================================================================
#
# Some defaults variables
#
# ========================================================================
plt.rc("text", usetex=True)
plt.rc("figure", max_open_warning=100)
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


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == "__main__":

    # Setup
    ct.add_directory(".")

    # Define fuel and oxidizer
    fuel = ct.Solution("LiDryer.cti", "gas")
    fuel_XH2 = 0.45
    fuel_XN2 = 0.55
    fuel.TPX = 300.0, ct.one_atm, f"H2:{fuel_XH2},N2:{fuel_XN2}"

    oxidizer = ct.Solution("LiDryer.cti", "gas")
    oxidizer_XO2 = 0.21
    oxidizer_XN2 = 0.79
    oxidizer.TPX = 300.0, ct.one_atm, f"O2:{oxidizer_XO2},N2:{oxidizer_XN2}"

    # Define blending
    T0 = 300.0
    Lx = 2 * np.pi * 1e-3
    R = 0.785 * 1e-3
    c = 3
    N = 2048
    x = np.linspace(0, Lx, N)
    Rd = np.sqrt((x - 0.5 * Lx) ** 2)
    phi = 0.5 * (1 + np.tanh(c * (Rd - R) / R))
    YH2s = fuel.Y[fuel.species_index("H2")] * (1 - phi)
    YO2s = oxidizer.Y[oxidizer.species_index("O2")] * phi

    # Loop on each cell
    pressure = np.zeros(len(x))
    temp = np.zeros(len(x))
    YH2s_equ = np.zeros(len(x))
    YO2s_equ = np.zeros(len(x))
    YH2Os_equ = np.zeros(len(x))
    YHs_equ = np.zeros(len(x))
    YOs_equ = np.zeros(len(x))
    YOHs_equ = np.zeros(len(x))
    YHO2s_equ = np.zeros(len(x))
    YH2O2s_equ = np.zeros(len(x))
    YN2s_equ = np.zeros(len(x))
    for i, (YH2, YO2) in enumerate(zip(YH2s, YO2s)):
        mixture = ct.Solution("LiDryer.cti", "gas")
        mixture.TPY = T0, ct.one_atm, f"H2:{YH2},O2:{YO2},N2:{1-(YH2+YO2)}"

        mixture.equilibrate("HP")

        pressure[i] = mixture.P
        temp[i] = mixture.T
        YH2s_equ[i] = mixture.Y[mixture.species_index("H2")]
        YO2s_equ[i] = mixture.Y[mixture.species_index("O2")]
        YH2Os_equ[i] = mixture.Y[mixture.species_index("H2O")]
        YHs_equ[i] = mixture.Y[mixture.species_index("H")]
        YOs_equ[i] = mixture.Y[mixture.species_index("O")]
        YOHs_equ[i] = mixture.Y[mixture.species_index("OH")]
        YHO2s_equ[i] = mixture.Y[mixture.species_index("HO2")]
        YH2O2s_equ[i] = mixture.Y[mixture.species_index("H2O2")]
        YN2s_equ[i] = mixture.Y[mixture.species_index("N2")]

    # Save data
    df = pd.DataFrame(
        {
            "x": x,
            "YH2": YH2s,
            "YO2": YO2s,
            "YN2": 1 - (YH2s + YO2s),
            "Tadiabatic": temp,
            "pressure": pressure,
            "YH2equ": YH2s_equ,
            "YO2equ": YO2s_equ,
            "YN2equ": YN2s_equ,
            "YH2Oequ": YH2Os_equ,
            "YHequ": YHs_equ,
            "YOequ": YOs_equ,
            "YOHequ": YOHs_equ,
            "YHO2equ": YHO2s_equ,
            "YH2O2equ": YH2O2s_equ,
        }
    )
    df.to_csv("profiles.dat", index=False)

    fig1, ax1 = plt.subplots()
    p = ax1.plot(x, temp, ls="-", lw=2, color=cmap[0])
    p[0].set_dashes(dashseq[0])
    ax1.set_xlabel(r"$x~[cm]$", fontsize=22)
    ax1.set_ylabel(r"$T~[K]$", fontsize=22)
    ax1.set_ylim([200, 2000])
    plt.setp(ax1.get_xmajorticklabels(), fontsize=16)
    plt.setp(ax1.get_ymajorticklabels(), fontsize=16)

    ax2 = ax1.twinx()
    p = ax2.plot(x, YH2s, ls="-", lw=2, color=cmap[1])
    p[0].set_dashes(dashseq[1])
    p = ax2.plot(x, YO2s, ls="-", lw=2, color=cmap[2])
    p[0].set_dashes(dashseq[2])
    ax2.set_ylabel(r"$Y~[-]$", fontsize=22)
    ax2.set_ylim([0, 0.4])
    plt.setp(ax2.get_xmajorticklabels(), fontsize=16)
    plt.setp(ax2.get_ymajorticklabels(), fontsize=16)

    fig2, ax3 = plt.subplots()
    p = ax3.plot(x, temp, ls="-", lw=2, color=cmap[0])
    p[0].set_dashes(dashseq[0])
    ax3.set_xlabel(r"$x~[cm]$", fontsize=22)
    ax3.set_ylabel(r"$T~[K]$", fontsize=22)
    ax3.set_ylim([200, 2000])
    plt.setp(ax3.get_xmajorticklabels(), fontsize=16)
    plt.setp(ax3.get_ymajorticklabels(), fontsize=16)

    ax4 = ax3.twinx()
    p = ax4.plot(x, YH2s_equ, ls="-", lw=2, color=cmap[1])
    p[0].set_dashes(dashseq[1])
    p = ax4.plot(x, YO2s_equ, ls="-", lw=2, color=cmap[2])
    p[0].set_dashes(dashseq[2])
    ax4.set_ylabel(r"$Y~[-]$", fontsize=22)
    ax4.set_ylim([0, 0.3])
    plt.setp(ax4.get_xmajorticklabels(), fontsize=16)
    plt.setp(ax4.get_ymajorticklabels(), fontsize=16)

    # Save plots
    fname = "plots.pdf"
    with PdfPages(fname) as pdf:
        plt.figure(fig1.number)
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure(fig2.number)
        plt.tight_layout()
        pdf.savefig(dpi=300)
