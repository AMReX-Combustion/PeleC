#!/usr/bin/env python
#
# Generate a table of the velocity fluctuations for the homogeneous
# isotropic turbulence case at a specific k0 (default to 4)
#
# Order of operations:
#   1. velocity fluctuations generated on a 512^3 grid in wavenumber space
#   2. Coefficients associated to wavenumbers that cannot be represented on
#      the desired grid are set to 0 (sharp wavenumber cutoff)
#   3. inverse Fourier transform of the velocity fluctuations (512^3 grid)
#   4. velocity fluctuations resampled on the desired grid (N^3)
#
# The velocity fluctuations are normalized by urms0 so to get the
# actual velocity fluctuations, one must multiply these velocities by
# the appropriate urms0.
#
#

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import sys
import time
from datetime import timedelta
import numpy as np
import scipy.interpolate as spi
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt


# ========================================================================
#
# Parse arguments
#
# ========================================================================
parser = argparse.ArgumentParser(
    description="Generate the velocity fluctuations for the HIT IC"
)
parser.add_argument(
    "-k0", help="Wave number containing highest energy", type=float, default=4.0
)
parser.add_argument("-N", help="Resolution", type=int, default=16)
parser.add_argument(
    "-s", "--seed", help="Random number generator seed", type=int, default=42
)
parser.add_argument(
    "-p", "--plot", help="Save a plot of the x-velocity", action="store_true"
)
args = parser.parse_args()

# ===============================================================================
#
# Some defaults variables
#
# ===============================================================================
plt.rc("text", usetex=True)
plt.rc("font", family="serif", serif="Times")
cmap_med = [
    "#F15A60",
    "#7AC36A",
    "#5A9BD4",
    "#FAA75B",
    "#9E67AB",
    "#CE7058",
    "#D77FB4",
    "#737373",
]
cmap = [
    "#EE2E2F",
    "#008C48",
    "#185AA9",
    "#F47D23",
    "#662C91",
    "#A21D21",
    "#B43894",
    "#010202",
]
dashseq = [
    (None, None),
    [10, 5],
    [10, 4, 3, 4],
    [3, 3],
    [10, 4, 3, 4, 3, 4],
    [3, 3],
    [3, 3],
]
markertype = ["s", "d", "o", "p", "h"]

# ========================================================================
#
# Function definitions
#
# ========================================================================
def div0(a, b):
    """Ignore division by 0, just replace it by 0,

    From: http://stackoverflow.com/questions/26248654/numpy-return-0-with-divide-by-zero
    e.g. div0( [-1, 0, 1], 0 ) -> [0, 0, 0]
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        c = np.true_divide(a, b)
        c[~np.isfinite(c)] = 0  # -inf inf NaN
    return c


# ========================================================================
def abs2(x):
    """This is equivalent to np.abs(x)**2 or x*np.conj(x)

    To make it faster, add this right before the function definition
    import numba
    @numba.vectorize([numba.float64(numba.complex128),numba.float32(numba.complex64)])
    """
    return x.real**2 + x.imag**2


# ========================================================================
#
# Main
#
# ========================================================================

# Timer
start = time.time()

# ========================================================================
# 1. velocity fluctuations generated on a 512^3 grid in wavenumber space

# Dimension of the large cube
N = 512
halfN = int(round(0.5 * N))
xs = 0
xe = 2.0 * np.pi
L = xe - xs
dx = L / N

# Only work if N and args.N are even
if not ((args.N % 2 == 0) and N % 2 == 0):
    print("N or args.N is not even. Exiting")
    sys.exit(1)

# Get cell centered values and meshed grid
x = np.linspace(xs, xe, N + 1)
xc = (x[1:] + x[:-1]) / 2  # get cell center coordinates
X, Y, Z = np.meshgrid(xc, xc, xc, indexing="ij")

# Get the wave numbers and associated quantities
k = np.concatenate((np.arange(halfN), np.arange(-halfN, 0, 1)), axis=0)
khalf = np.arange(halfN + 1)
k1, k2, k3 = np.meshgrid(k, k, khalf, indexing="ij")
kmag = np.sqrt(k1**2 + k2**2 + k3**2)
k12 = np.sqrt(k1**2 + k2**2)
k1k12 = div0(k1, k12)
k2k12 = div0(k2, k12)
k3kmag = div0(k3, kmag)
k12kmag = div0(k12, kmag)

# Generate data

# # Toy Fourier data corresponding to uo = cos(X) * cos(2*Y) * cos(3*Z)
# uo = np.cos(X) * np.cos(2*Y) * np.cos(3*Z)
# uf = np.fft.rfftn(uo)
# vf = np.copy(uf)
# wf = np.copy(uf)

# Energy spectrum
Ek = (
    16.0
    * np.sqrt(2.0 / np.pi)
    * (kmag**4)
    / (args.k0**5)
    * np.exp(-2.0 * (kmag**2) / (args.k0**2))
)

# Draw random numbers
np.random.seed(args.seed)
phi1 = np.random.uniform(0, 2 * np.pi, np.shape(kmag))
phi2 = np.random.uniform(0, 2 * np.pi, np.shape(kmag))
phi3 = np.random.uniform(0, 2 * np.pi, np.shape(kmag))

# the random quantities
prefix = np.sqrt(2.0 * div0(Ek, 4.0 * np.pi * (kmag**2)))
a = prefix * np.exp(1j * phi1) * np.cos(phi3)
b = prefix * np.exp(1j * phi2) * np.sin(phi3)

# the random velocities
uf = k2k12 * a + k1k12 * k3kmag * b
vf = k2k12 * k3kmag * b - k1k12 * a
wf = -k12kmag * b

# Impose the 3D spherical symmetry (to ensure we have a real signal)
# equiv: uf[-l,-m,0] = np.conj(uf[ l, m,0]) for l=0..N/2 and m=0..N/2
uf[N:halfN:-1, N:halfN:-1, 0] = np.conj(uf[1:halfN, 1:halfN, 0])
# symmetry on first column
uf[N:halfN:-1, 0, 0] = np.conj(uf[1:halfN, 0, 0])
# symmetry on first row
uf[0, N:halfN:-1, 0] = np.conj(uf[0, 1:halfN, 0])
# symmetry about the (halfN,halfN) element
uf[halfN - 1 : 0 : -1, N : halfN - 1 : -1, 0] = np.conj(
    uf[halfN + 1 : N, 1 : halfN + 1, 0]
)

vf[N:halfN:-1, N:halfN:-1, 0] = np.conj(vf[1:halfN, 1:halfN, 0])
vf[halfN - 1 : 0 : -1, N : halfN - 1 : -1, 0] = np.conj(
    vf[halfN + 1 : N, 1 : halfN + 1, 0]
)
vf[N:halfN:-1, 0, 0] = np.conj(vf[1:halfN, 0, 0])
vf[0, N:halfN:-1, 0] = np.conj(vf[0, 1:halfN, 0])

wf[N:halfN:-1, N:halfN:-1, 0] = np.conj(wf[1:halfN, 1:halfN, 0])
wf[halfN - 1 : 0 : -1, N : halfN - 1 : -1, 0] = np.conj(
    wf[halfN + 1 : N, 1 : halfN + 1, 0]
)
wf[N:halfN:-1, 0, 0] = np.conj(wf[1:halfN, 0, 0])
wf[0, N:halfN:-1, 0] = np.conj(wf[0, 1:halfN, 0])

# Normalize. Because we are generating the data in wavenumber space,
# we have to multiply by N**3 because in the definition of the numpy
# ifftn there is a 1/N**n.
uf = uf * N**3
vf = vf * N**3
wf = wf * N**3

# # Quick check on energy content (make sure you add both the current
# # contribution and the one we are neglecting because we are assuming
# # real input data)
# print('Energy = int E(k) dk = 0.5 * int (uf**2 + vf**2 wf**2) dk1 dk2 dk3 = {0:.10f} ~= 3/2'.format(
#     (np.sum(abs2(uf          ) + abs2(vf          ) + abs2(wf          )) +
#      np.sum(abs2(uf[:,:,1:-1]) + abs2(vf[:,:,1:-1]) + abs2(wf[:,:,1:-1])))
#     * 0.5 / N**6))

# if plotting, save the original field (before filtering)
if args.plot:
    uo = np.fft.irfftn(uf)
    Eko = (
        16.0
        * np.sqrt(2.0 / np.pi)
        * (khalf**4)
        / (args.k0**5)
        * np.exp(-2.0 * (khalf**2) / (args.k0**2))
    )

    # Get the spectrum from 3D velocity field
    kbins = np.arange(1, halfN + 1)
    Nbins = len(kbins)
    whichbin = np.digitize(kmag.flat, kbins)
    ncount = np.bincount(whichbin)

    KI = (abs2(uf) + abs2(vf) + abs2(wf)) * 0.5 / N**6
    KI[:, :, 1:-1] += (
        (abs2(uf[:, :, 1:-1]) + abs2(vf[:, :, 1:-1]) + abs2(wf[:, :, 1:-1]))
        * 0.5
        / N**6
    )

    Eku = np.zeros(len(ncount) - 1)
    for n in range(1, len(ncount)):
        Eku[n - 1] = np.sum(KI.flat[whichbin == n])

    ku = 0.5 * (kbins[0 : Nbins - 1] + kbins[1:Nbins]) + 1
    Eku = Eku[1:Nbins]


# ========================================================================
# 2. Coefficients associated to wavenumbers that cannot be represented
# on the desired grid are set to 0 (sharp wavenumber cutoff)
kmagc = 0.5 * args.N
uf[kmag > kmagc] = 0.0
vf[kmag > kmagc] = 0.0
wf[kmag > kmagc] = 0.0


# ========================================================================
# 3. inverse Fourier transform of the velocity fluctuations (512^3 grid)
u = np.fft.irfftn(uf, s=(N, N, N))
v = np.fft.irfftn(vf, s=(N, N, N))
w = np.fft.irfftn(wf, s=(N, N, N))

# Another energy content check
print(
    "Energy = 1/V * int E(x,y,z) dV = 0.5/V * int (u**2 + v**2 + w**2) dx dy dz = {0:.10f} ~= 3/2".format(
        np.sum(u**2 + v**2 + w**2) * 0.5 * (dx / L) ** 3
    )
)

# # Enstrophy check
# _, dudy, dudz = np.gradient(u, dx)
# dvdx, _, dvdz = np.gradient(v, dx)
# dwdx, dwdy, _ = np.gradient(w, dx)
# wx = dwdy-dvdz
# wy = dudz-dwdx
# wz = dvdx-dudy
# lambda0 = 2.0/args.k0
# print('Enstrophy = 0.5/V * int (wx**2 + wy**2 + wz**2) dx dy dz=
# {0:.10f} ~= '.format(np.sum(wx**2+wy**2+wz**2) * 0.5 * (dx/L)**3 *
# lambda0**2))

# ========================================================================
# 4. velocity fluctuations re-sampled on the desired grid (N^3)
xr = np.linspace(xs, xe, args.N + 1)
xrc = (xr[1:] + xr[:-1]) / 2
Xr, Yr, Zr = np.meshgrid(xrc, xrc, xrc, indexing="ij")

Xr = Xr.reshape(-1, order="F")
Yr = Yr.reshape(-1, order="F")
Zr = Zr.reshape(-1, order="F")

ur = spi.interpn((xc, xc, xc), u, (Xr, Yr, Zr), method="linear")
vr = spi.interpn((xc, xc, xc), v, (Xr, Yr, Zr), method="linear")
wr = spi.interpn((xc, xc, xc), w, (Xr, Yr, Zr), method="linear")


# ========================================================================
# Save the data in Fortran ordering
fname = "hit_ic_{0:d}_{1:d}.dat".format(int(args.k0), args.N)
data = np.vstack((Xr, Yr, Zr, ur, vr, wr)).T
np.savetxt(fname, data, fmt="%.18e", delimiter=",", header="x, y, z, u, v, w")


# ========================================================================
# plot (only u fluctuations)
if args.plot:
    import matplotlib as mpl

    mpl.use("Agg")
    import matplotlib.pyplot as plt

    datmin = u.min()
    datmax = u.max()
    # print("min/max u =",datmin,datmax)

    # Original data
    # transpose and origin change bc I used meshgrid with ij and not xy
    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(14, 14))
    ax[0, 0].imshow(
        uo[:, :, 0].T,
        origin="lower",
        extent=[xs, xe, xs, xe],
        cmap="RdBu_r",
        vmin=datmin,
        vmax=datmax,
    )
    ax[0, 0].set_title("Original data (x,y)")
    ax[0, 1].imshow(
        uo[:, 0, :].T,
        origin="lower",
        extent=[xs, xe, xs, xe],
        cmap="RdBu_r",
        vmin=datmin,
        vmax=datmax,
    )
    ax[0, 1].set_title("Original data (x,z)")
    ax[0, 2].imshow(
        uo[0, :, :].T,
        origin="lower",
        extent=[xs, xe, xs, xe],
        cmap="RdBu_r",
        vmin=datmin,
        vmax=datmax,
    )
    ax[0, 2].set_title("Original data (y,z)")

    # Filtered original data
    ax[1, 0].imshow(
        u[:, :, 0].T,
        origin="lower",
        extent=[xs, xe, xs, xe],
        cmap="RdBu_r",
        vmin=datmin,
        vmax=datmax,
    )
    ax[1, 0].set_title("Filtered original data (x,y)")
    ax[1, 1].imshow(
        u[:, 0, :].T,
        origin="lower",
        extent=[xs, xe, xs, xe],
        cmap="RdBu_r",
        vmin=datmin,
        vmax=datmax,
    )
    ax[1, 1].set_title("Filtered original data (x,z)")
    ax[1, 2].imshow(
        u[0, :, :].T,
        origin="lower",
        extent=[xs, xe, xs, xe],
        cmap="RdBu_r",
        vmin=datmin,
        vmax=datmax,
    )
    ax[1, 2].set_title("Filtered original data (y,z)")

    # Downsampled filtered data
    ur = ur.reshape(args.N, args.N, args.N, order="F")
    ax[2, 0].imshow(
        ur[:, :, 0].T,
        origin="lower",
        extent=[xs, xe, xs, xe],
        cmap="RdBu_r",
        vmin=datmin,
        vmax=datmax,
    )
    ax[2, 0].set_title("Downsampled data (x,y)")
    ax[2, 1].imshow(
        ur[:, 0, :].T,
        origin="lower",
        extent=[xs, xe, xs, xe],
        cmap="RdBu_r",
        vmin=datmin,
        vmax=datmax,
    )
    ax[2, 1].set_title("Downsampled data (x,z)")
    ax[2, 2].imshow(
        ur[0, :, :].T,
        origin="lower",
        extent=[xs, xe, xs, xe],
        cmap="RdBu_r",
        vmin=datmin,
        vmax=datmax,
    )
    ax[2, 2].set_title("Downsampled data (y,z)")

    plt.savefig("hit_ic_u_{0:d}_{1:d}.png".format(int(args.k0), args.N), format="png")

    # Fourier coefficients of original data
    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(14, 8))
    ax[0, 0].imshow(np.real(uf[:, :, 0].T), origin="lower", cmap="RdBu_r")
    ax[0, 0].set_title("Real Fourier coefficients (x,y)")
    ax[0, 1].imshow(np.real(uf[:, 0, :].T), origin="lower", cmap="RdBu_r")
    ax[0, 1].set_title("Real Fourier coefficients (x,z)")
    ax[0, 2].imshow(np.real(uf[0, :, :].T), origin="lower", cmap="RdBu_r")
    ax[0, 2].set_title("Real Fourier coefficients (y,z)")
    ax[1, 0].imshow(np.imag(uf[:, :, 0].T), origin="lower", cmap="RdBu_r")
    ax[1, 0].set_title("Imag Fourier coefficients (x,y)")
    ax[1, 1].imshow(np.imag(uf[:, 0, :].T), origin="lower", cmap="RdBu_r")
    ax[1, 1].set_title("Imag Fourier coefficients (x,z)")
    ax[1, 2].imshow(np.imag(uf[0, :, :].T), origin="lower", cmap="RdBu_r")
    ax[1, 2].set_title("Imag Fourier coefficients (y,z)")
    plt.savefig("hit_ic_uf_{0:d}_{1:d}.png".format(int(args.k0), args.N), format="png")

    # Spectrum
    plt.figure(20)
    ax = plt.gca()
    p = plt.loglog(khalf, Eko, color=cmap[-1], lw=2)
    p[0].set_dashes(dashseq[0])
    p = plt.loglog(ku, Eku, color=cmap[0], lw=2)
    p[0].set_dashes(dashseq[1])
    plt.ylim([1e-16, 10])
    plt.xlabel(r"$k$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$E(k)$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    plt.savefig(
        "hit_ic_spectrum_{0:d}_{1:d}.png".format(int(args.k0), args.N), format="png"
    )

# output timer
end = time.time() - start
print("Elapsed time " + str(timedelta(seconds=end)) + " (or {0:f} seconds)".format(end))
