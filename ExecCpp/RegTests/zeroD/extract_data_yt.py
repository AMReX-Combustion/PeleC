import yt
from sys import argv
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import os

# ==================================================
# this script extracts average temperature
# from multiple plot files
# example use python extract_data_yt.py "plt?????"
# ==================================================

fn_pattern = argv[1]
fn_list = sorted(glob.glob(fn_pattern))

print(fn_list)
fname = "zerod_soln.dat"
outfile = open(fname, "w")

for i, fn in enumerate(fn_list):
    ds = yt.load(fn)
    ad = ds.all_data()
    meanTemp = np.mean(ad["Temp"])
    outfile.write("%e\t%e\n" % (ds.current_time, meanTemp))

outfile.close()

# plot
ref = pd.read_csv("cantera_soln", delim_whitespace=True)
df = pd.read_csv(fname, delim_whitespace=True, header=None, names=["t", "temp"])

plt.figure()
plt.plot(ref["Time(s)"], ref["Temperature(K)"], label="Ref.")
plt.plot(df["t"], df["temp"], label="PeleC")
plt.legend()
plt.savefig("plot.png")
