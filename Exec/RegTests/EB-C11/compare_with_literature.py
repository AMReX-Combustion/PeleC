import yt
from sys import argv
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# =======================================
# read literature solutions
def readlitsoln(fname):

    infile = open(fname, "r")
    x = np.array([])
    d = np.array([])

    for line in infile:
        x = np.append(x, float(line.split()[0]))
        d = np.append(d, float(line.split()[1]))

    return (x, d)


# =======================================


ds = yt.load(argv[1])

L = (ds.domain_right_edge - ds.domain_left_edge).d
axial_dir = np.argmax(L)

clength = L[axial_dir]
cdepth = L[(axial_dir + 1) % 3]

fieldname = "density"
res = 100
slicedepth = cdepth / 2
slc = ds.slice((axial_dir + 2) % 3, slicedepth)
frb = slc.to_frb((clength, "cm"), res)
x = np.linspace(0, clength, res)
fld = np.array(frb[fieldname])[int(res / 2), :]

outfile = open("pelec_soln.dat", "w")
for i in range(len(x)):
    outfile.write("%e\t%e\n" % (x[i], fld[i]))
outfile.close()

fname1 = "Lv_Ihme_JCP_2014"
fname2 = "Grogan_Ihme_ShockWaves_2020"

(x1, d1) = readlitsoln(fname1)
(x2, d2) = readlitsoln(fname2)
# =======================================

print("plotting")

plt.plot(x, fld / np.max(fld), "k", label="PeleC")
plt.plot(x1, d1, "r*-", label=fname1 + " (N2 left)")
plt.plot(x2, d2, "g*-", label=fname2 + " (HE left)")
plt.legend(loc="best")
plt.suptitle("Multi-species shock tube density solution")

plt.savefig("dens_shock_tube.png")
# =======================================
