import yt
from sys import argv
import numpy as np
import glob
import pandas as pd

# ==================================================
# this script extracts average temperature
# from multiple plot files
# example use python extract_data_yt.py "plt?????"
# ==================================================

fn_pattern = argv[1]
outfile = argv[2]
fn_list = sorted(glob.glob(fn_pattern))

print(fn_list)
times = []
minTemps = []
avgTemps = []
maxTemps = []

for i, fn in enumerate(fn_list):
    ds = yt.load(fn)
    ad = ds.all_data()
    Temp = np.array(ad["Temp"])[np.where(ad["vfrac"] > 1e-8)]
    minTemp = np.min(Temp)
    avgTemp = np.sum(ad["Temp"] * ad["vfrac"]) / np.sum(ad["vfrac"])
    maxTemp = np.max(Temp)

    times.append(float(ds.current_time))
    minTemps.append(minTemp)
    avgTemps.append(float(avgTemp))
    maxTemps.append(maxTemp)

    print(fn, float(ds.current_time), minTemp, avgTemp, maxTemp)

dataframe = pd.DataFrame(
    {
        "file": fn_list,
        "time": times,
        "min T": minTemps,
        "max T": maxTemps,
        "avg T": avgTemps,
    }
)

print(dataframe)
dataframe.to_csv(outfile)
