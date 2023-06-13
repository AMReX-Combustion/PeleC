import os
import glob
import argparse
import yt
import pandas as pd
import numpy as np


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

    fields = [
        "x",
        "y",
        "z",
        "density",
        "pressure",
        "x_velocity",
        "y_velocity",
        "z_velocity",
        "Temp",
        "vfrac",
    ]
    angle = 30.0
    length = 4.0
    factor = 0.05
    tube_half_height = length * factor
    cost = np.cos(np.radians(angle))
    sint = np.sin(np.radians(angle))
    tant = np.tan(np.radians(angle))
    for fdir in args.fdirs:
        pltfiles = sorted(glob.glob(os.path.join(fdir, "plt*")))

        ds = yt.load(pltfiles[-1])

        start = [0.0, 0.0, 0.0]
        end = [length, length * tant, 0.0]
        ray = ds.r[start:end]
        srt = np.argsort(ray[("boxlib", "x")])

        df = pd.DataFrame({f: np.array(ray[("boxlib", f)][srt]) for f in fields})
        df.Temp /= 3.48429e-07
        df["time"] = ds.current_time
        df["levels"] = ds.max_level + 1
        df["dxmin"] = ds.index.get_smallest_dx()
        df["xp"] = df.x * cost + df.y * sint - 0.5 * length
        df["yp"] = df.x * sint - df.y * cost
        df["u"] = df.x_velocity * cost + df.y_velocity * sint
        df["v"] = -df.x_velocity * sint + df.y_velocity * cost
        df.sort_values(by=["xp"], inplace=True)
        df.to_csv(os.path.join(fdir, "profiles.csv"), index=False)
