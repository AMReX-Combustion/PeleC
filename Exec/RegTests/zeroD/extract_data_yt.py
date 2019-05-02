import yt
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

#==================================================
#this script extracts average temperature
#from multiple plot files
#example use python extract_data_yt.py "plt?????"
#==================================================

fn_pattern = argv[1]
fn_list = sorted(glob.glob(fn_pattern), key=os.path.getmtime)

print(fn_list)
outfile=open("zerod_soln.dat","w")

for i, fn in enumerate(fn_list):
    ds=yt.load(fn)
    ad=ds.all_data()
    meanTemp=np.mean(ad["Temp"])
    outfile.write("%e\t%e\n"%(ds.current_time,meanTemp))

outfile.close()
