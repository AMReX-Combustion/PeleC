import os
import glob
import argparse
import yt
import pandas as pd
import numpy as np


def read_ke(fname):
    lst = []
    with open(fname, "r") as f:
        for line in f:
            if line.startswith("TIME= ") and "RHO*K" in line:
                line = line.split()
                lst.append([float(line[1]), float(line[-1])])

    return pd.DataFrame(np.array(lst), columns=["time", "rho_ke"])


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

    fields = ["y", "x_velocity", "y_velocity", "z_velocity", "vfrac"]
    xlocs = [0, 3, 6, 9, 12 - 1e-10]
    for fdir in args.fdirs:
        pltfiles = sorted(glob.glob(os.path.join(fdir, "plt*")))

        # Get time history of KE
        df = read_ke(os.path.join(fdir, "out"))
        df.to_csv(os.path.join(fdir, "history.csv"), index=False)

        # Get the profiles
        ds = yt.load(pltfiles[-1])
        lst = []
        for xloc in xlocs:
            ray = ds.ortho_ray(1, (0, xloc))
            srt = np.argsort(ray[("boxlib", "y")])

            df = pd.DataFrame({f: np.array(ray[("boxlib", f)][srt]) for f in fields})
            df["xloc"] = xloc
            df = df[df.vfrac > 0]
            lst.append(df)

        df = pd.concat(lst)
        df["time"] = ds.current_time
        df.to_csv(os.path.join(fdir, "profiles.csv"), index=False)
