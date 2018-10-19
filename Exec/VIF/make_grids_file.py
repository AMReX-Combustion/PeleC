#!/usr/bin/env python3
#
# This file creates a grids file for the VIF case
#

# ========================================================================
#
# Imports
#
# ========================================================================
import numpy as np


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == '__main__':

    # ========================================================================
    # Setup

    # ratio of Ly to Lx
    L = 10

    # number of cells on the coarse level
    nyc = 64
    nxc = nyc * L

    # Start and end in x direction of fine zone
    Lxs = int(0 * nyc / 2)
    Lxe = int(12 * nyc / 2)

    # Grid refinement ratio
    ratio = 2

    # number of levels (only 1 is supported right now)
    nlvls = 1

    # max grid size on each level
    max_grid_size = 64

    # output file name
    oname = 'grids_file_{0:d}'.format(nyc)

    # ========================================================================
    # Write output file
    with open(oname, 'w') as of:
        of.write(str(nlvls) + '\n')

        xs = np.arange(Lxs,
                       Lxe,
                       int(max_grid_size))
        xe = xs + int(max_grid_size) - 1

        ys = np.arange(0,
                       int(nyc),
                       int(max_grid_size))

        ye = ys + int(max_grid_size) - 1

        XS = np.meshgrid(xs, ys)
        XE = np.meshgrid(xe, ye)
        of.write('{0:d}\n'.format(len(XS[0].flatten())))

        for xs, ys, xe, ye in zip(XS[0].flatten(),
                                  XS[1].flatten(),
                                  XE[0].flatten(),
                                  XE[1].flatten()):
            of.write("(({0:d},{1:d})({2:d},{3:d}))\n".format(
                xs, ys, xe, ye))

        # need a blank line at the end of the grids file
        of.write('\n')
