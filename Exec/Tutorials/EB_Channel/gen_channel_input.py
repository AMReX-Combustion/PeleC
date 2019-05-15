#!/usr/bin/env python3

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import math


# ========================================================================
#
# Functions
#
# ========================================================================
def round_up_to_multiple(num, divisor):
    return num - (num % divisor) + divisor


# ========================================================================
def write_probin(pname, reynolds, mach, prandtl, angle):
    with open(pname, "w") as f:

        f.write("&fortin\n")
        f.write("\n")
        f.write(f"  reynolds = {reynolds}\n")
        f.write(f"  mach = {mach}\n")
        f.write(f"  prandtl = {prandtl}\n")
        f.write(f"  angle = {angle}\n")
        f.write("\n")
        f.write("/\n")
        f.write("\n")
        f.write("&tagging\n")
        f.write("\n")
        f.write("  denerr = 1e20\n")
        f.write("  dengrad = 1e20\n")
        f.write("  max_denerr_lev = 3\n")
        f.write("  max_dengrad_lev = 3\n")
        f.write("\n")
        f.write("  presserr = 1e20\n")
        f.write("  pressgrad = 1e20\n")
        f.write("  max_presserr_lev = 3\n")
        f.write("  max_pressgrad_lev = 3\n")
        f.write("  \n")
        f.write("  velerr = 1e20\n")
        f.write("  velgrad = 1e20\n")
        f.write("  max_velerr_lev = 3\n")
        f.write("  max_velgrad_lev = 3\n")
        f.write("\n")
        f.write("  temperr  = 1e20\n")
        f.write("  tempgrad = 1e20\n")
        f.write("  max_velerr_lev = 5\n")
        f.write("  max_velgrad_lev = 5\n")
        f.write("\n")
        f.write("  vorterr = 1000\n")
        f.write("  max_vorterr_lev = 5\n")
        f.write("  \n")
        f.write("  vfracerr = 1e20\n")
        f.write("  max_vfracerr_lev = 5\n")
        f.write("\n")
        f.write("/\n")
        f.write("\n")
        f.write("&extern\n")
        f.write("  eos_gamma = 1.4\n")
        f.write("\n")
        f.write("\n")
        f.write("/\n")


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(
        description="A simple tool to the channel flow setup"
    )
    parser.add_argument(
        "-a", "--angle", help="Flow angle (degrees)", type=float, required=True
    )
    parser.add_argument(
        "-r", "--reynolds", help="Reynolds number", type=float, required=True
    )
    parser.add_argument("--mach", help="Mach number", type=float, default=0.1)
    parser.add_argument("--prandtl", help="Prandtl number", type=float, default=0.71)
    parser.add_argument(
        "-n",
        "--ncell",
        help="Number of cells per channel half height",
        type=float,
        required=True,
        choices=[8 * i for i in range(1, 1001)],
    )
    parser.add_argument(
        "-t",
        "--translation",
        help="Translation of channel in y-direction",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--factor", help="Multiplier for the laminar length", type=float, default=2.0
    )
    args = parser.parse_args()

    # Laminar channel flow development
    channel_half_height = 1.0
    laminar_length = 0.05 * args.reynolds * 2 * channel_half_height
    domain_length = args.factor * laminar_length

    # Determine triangle coordinates
    m = math.tan(math.radians(args.angle))
    p = channel_half_height / math.cos(math.radians(args.angle)) + args.translation
    buf = 10 * domain_length

    def bot(x):
        return m * x - p

    def top(x):
        return m * x + p

    p0b = [-buf, bot(-buf)]
    p1b = [domain_length + buf, -buf]
    p2b = [domain_length + buf, bot(domain_length + buf)]

    p0t = [-buf, top(-buf)]
    p1t = [domain_length + buf, top(domain_length + buf)]
    p2t = [-buf, buf]

    tri_line = (
        f"extruded_triangles.tri_0_point_0={p0b[0]} {p0b[1]} 0.0 "
        f"extruded_triangles.tri_0_point_1={p1b[0]} {p1b[1]} 0.0 "
        f"extruded_triangles.tri_0_point_2={p2b[0]} {p2b[1]} 0.0 "
        f"extruded_triangles.tri_1_point_0={p0t[0]} {p0t[1]} 0.0 "
        f"extruded_triangles.tri_1_point_1={p1t[0]} {p1t[1]} 0.0 "
        f"extruded_triangles.tri_1_point_2={p2t[0]} {p2t[1]} 0.0 "
    )

    # Problem geometry
    dx = channel_half_height / args.ncell
    blocking_factor = 8
    nbyp = int(round_up_to_multiple(1.05 * top(domain_length) / dx, blocking_factor))
    nbym = int(blocking_factor * 2)
    total_ncelly = int(2 * args.ncell + nbyp + nbym)
    total_ncellx = int(domain_length / dx)
    total_ncellz = blocking_factor

    lo = [0.0, -(nbym + args.ncell) * dx, 0.0]
    hi = [total_ncellx * dx, (nbyp + args.ncell) * dx, total_ncellz * dx]

    geom = (
        f"geometry.prob_lo={lo[0]} {lo[1]} {lo[2]} "
        f"geometry.prob_hi={hi[0]} {hi[1]} {hi[2]} "
        f"amr.n_cell={total_ncellx} {total_ncelly} {total_ncellz} "
    )

    # Generate probin file
    pname = "probin.channel"
    write_probin(pname, args.reynolds, args.mach, args.prandtl, args.angle)
    probin_line = f"amr.probin_file={pname}"

    # Construct the input line
    print(geom, tri_line, probin_line)
