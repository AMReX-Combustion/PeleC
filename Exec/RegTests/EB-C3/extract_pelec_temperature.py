import yt
from sys import argv
import numpy as np
import glob

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
    meanTemp = np.mean(ad["Temp"] * ad["vfrac"])
    meanvfrac = np.mean(ad["vfrac"])
    outfile.write("%e\t%e\n" % (ds.current_time, meanTemp / meanvfrac))

outfile.close()
