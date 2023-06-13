import yt
from sys import argv
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ds = yt.load(argv[1])
axial_dir = int(argv[2])

clength = 1.0
cwidth = 0.0625
cdepth = 0.0625

fieldname = "density"
res = 100
slicedepth = cdepth / 2
slc = ds.slice((axial_dir + 2) % 3, slicedepth)
frb = slc.to_frb((1, "cm"), res)
x = np.linspace(0, 1, res)
fld = np.array(frb[fieldname])[int(res / 2), :]

outfile = open("pelec_soln.dat", "w")

for i in range(len(x)):
    outfile.write("%e\t%e\n" % (x[i], fld[i]))

outfile.close()
