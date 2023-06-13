import os
import glob
import argparse
import yt
import pandas as pd
import numpy as np
import re


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
    return sorted(l, key=alphanum_key)


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

    fields = ["x", "density", "x_velocity", "pressure"]
    idxs = [0, -1]
    for fdir in args.fdirs:
        pltfiles = natural_sort(glob.glob(os.path.join(fdir, "plt*")))

        lst = []
        for idx in idxs:
            ds = yt.load(pltfiles[idx])
            ray = ds.ortho_ray(0, (0, 0))
            srt = np.argsort(ray[("boxlib", "x")])

            df = pd.DataFrame({f: np.array(ray[("boxlib", f)][srt]) for f in fields})
            df["time"] = ds.current_time
            lst.append(df)

        df = pd.concat(lst)
        df.to_csv(os.path.join(fdir, "profiles.csv"), index=False)
