import yt
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
    
    
fn_pattern = argv[1]
fn_list = sorted(glob.glob(fn_pattern), key=lambda f: int(f.split("plt")[1]))

print(fn_list)

minmaxset=False
if(len(argv) > 3):
    minval=float(argv[3])
    maxval=float(argv[4])
    minmaxset=True

for i, fn in enumerate(fn_list):

    ds=yt.load(fn)
    ad=ds.smoothed_covering_grid(level=ds.index.max_level,left_edge=[0.0,0.0,0.0],\
            dims=ds.domain_dimensions * 2**(ds.index.max_level))
    dxmin = ds.index.get_smallest_dx()

    dens = np.array(ad["density"])
    vfrac = np.array(ad["vfrac"])

    nz=np.shape(dens)[2]

    dengradvec=np.gradient(dens[:,:,int(nz/2)],dxmin)
    dengrad=np.sqrt(dengradvec[0]**2+dengradvec[1]**2)
    np.shape(dengrad)

    vfrac2d=vfrac[:,:,int(nz/2)]

    dengrad[vfrac2d==0.0]=np.nan

    fig0 = plt.figure(0)
    ax0 = fig0.add_subplot(111)
    im = ax0.imshow(
            dengrad.T,
            origin="lower",
            cmap="gray",
            aspect="equal"
            ,vmin=0.0,vmax=0.002
        )
    cbar = plt.colorbar(im, ax=ax0,orientation="horizontal")
    #ax0.contour(dengrad[0], levels=np.array([Zst]), colors='white', alpha=0.5)
    ax0.axis('off')
    cbar.ax.set_title("density gradient")
    plt.savefig("dengrad_%3.3d.png"%(i))
