#!/usr/bin/env python3

import argparse
import math


def round_up_to_multiple(num, divisor):
    return num - (num % divisor) + divisor


def get_prob(state, angle):
    return f"prob.pl={state[0]} prob.rhol={state[1]} prob.pr={state[2]} prob.rhor={state[3]} prob.angle={angle}"


if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(
        description="A simple tool to the shocktube flow setup"
    )
    parser.add_argument(
        "-a", "--angle", help="Flow angle (degrees)", type=float, required=True
    )
    parser.add_argument(
        "-s",
        "--state",
        nargs="+",
        help="Left-right state [pl, rhol, pr, rhor]",
        type=float,
        default=[1.0, 1.0, 0.1, 0.125],
    )
    parser.add_argument(
        "-n",
        "--ncell",
        help="Number of cells per tube half height",
        type=float,
        required=True,
        choices=[8 * i for i in range(1, 1001)],
    )
    parser.add_argument(
        "-t",
        "--translation",
        help="Translation of tube in y-direction",
        type=float,
        default=0.0,
    )
    parser.add_argument("--length", help="Domain length", type=float, default=4.0)
    args = parser.parse_args()

    # Tube description
    args.length = args.length
    factor = 0.05
    tube_half_height = factor * args.length

    # Determine triangle coordinates
    m = math.tan(math.radians(args.angle))
    p = tube_half_height / math.cos(math.radians(args.angle)) + args.translation
    buf = 10 * args.length

    def bot(x):
        return m * x - p

    def top(x):
        return m * x + p

    p0b = [-buf, bot(-buf)]
    p1b = [args.length + buf, -buf]
    p2b = [args.length + buf, bot(args.length + buf)]

    p0t = [-buf, top(-buf)]
    p1t = [args.length + buf, top(args.length + buf)]
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
    dx = tube_half_height / args.ncell
    blocking_factor = 8
    nbyp = int(round_up_to_multiple(1.05 * top(args.length) / dx, blocking_factor))
    nbym = int(blocking_factor * 2)
    total_ncelly = int(2 * args.ncell + nbyp + nbym)
    total_ncellx = int(args.length / dx)
    total_ncellz = blocking_factor

    lo = [0.0, -(nbym + args.ncell) * dx, 0.0]
    hi = [total_ncellx * dx, (nbyp + args.ncell) * dx, total_ncellz * dx]

    geom = (
        f"geometry.prob_lo={lo[0]} {lo[1]} {lo[2]} "
        f"geometry.prob_hi={hi[0]} {hi[1]} {hi[2]} "
        f"amr.n_cell={total_ncellx} {total_ncelly} {total_ncellz} "
    )

    # Generate probin file
    prob_line = get_prob(args.state, args.angle)

    # Construct the input line
    print(geom, tri_line, prob_line)
