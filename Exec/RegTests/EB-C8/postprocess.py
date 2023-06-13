import os
import glob
import argparse
import yt
import pandas as pd
import numpy as np


if __name__ == "__main__":
    fields = [
        "x",
        "y",
        "z",
        "density",
        "vfrac",
    ]
    angle = 30.0
    length = 0.87
    factor = 0.05
    tube_half_height = length * factor
    cost = np.cos(np.radians(angle))
    sint = np.sin(np.radians(angle))
    tant = np.tan(np.radians(angle))
    axislength = length / cost
    print("axis length:", axislength)

    outfile = open("densprofile.dat", "w")
    pltfiles = sorted(glob.glob(os.path.join(".", "plt*")))

    ds = yt.load(pltfiles[-1])

    start = [0.0, 0.0, 0.0]
    end = [length, length * tant, 0.0]
    ray = ds.r[start:end]
    srt = np.argsort(ray[("boxlib", "x")])

    df = pd.DataFrame({f: np.array(ray[("boxlib", f)][srt]) for f in fields})
    df["xp"] = (df.x * cost + df.y * sint) / axislength
    df["yp"] = df.x * sint - df.y * cost
    df.sort_values(by=["xp"], inplace=True)
    df.to_csv(
        os.path.join(".", "pelec_soln.dat"),
        sep="\t",
        columns=["xp", "density"],
        index=False,
    )
