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

    fields = ["x", "density"]
    idxs = [0, -1]
    for fdir in args.fdirs:
        pltfiles = sorted(glob.glob(os.path.join(fdir, "plt*")))

        lst = []
        for idx in idxs:
            ds = yt.load(pltfiles[idx])
            ray = ds.ortho_ray(0, (0, 0))
            srt = np.argsort(ray["x"])

            df = pd.DataFrame({f: np.array(ray[f][srt]) for f in fields})
            df["time"] = ds.current_time
            lst.append(df)

        df = pd.concat(lst)
        df.to_csv(os.path.join(fdir, "profiles.csv"), index=False)
