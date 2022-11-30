#
# Post-process data using VisIt. Generates data files of
# - the spatial average of u^+v^2+w^2
# - the spatial average of sound speed
# - the spatial average of u^2
# - the spatial average of Sij^2
# - peak temperature in time
# - lineouts of variables
#
# Usage:
#    visit -nowin -cli -s /path/to/pp_aux_vars.py -i inputs_3d
#

# ========================================================================
#
# Imports
#
# ========================================================================
import os
import sys
import re
import math
import time
from datetime import timedelta
import subprocess as sp
import argparse


# ========================================================================
#
# Function definitions
#
# ========================================================================
def truncate(number, digits):
    stepper = pow(10.0, digits)
    return math.trunc(stepper * number) / stepper


# ========================================================================
def get_simulation_params(icname, iname):
    with open(icname, "r") as f:
        line = [float(x) for x in f.read().split("\n")[1].split(",")]
        Lx = line[0]
        Ly = Lx

    with open(iname, "r") as f:
        for line in f:
            if "amr.n_cell" in line:
                ncells = [int(x) for x in re.findall(r"[0-9]+", line)]

    if len(ncells) == 3:
        Lz = Lx
    else:
        Lz = 0.0
    return Lx, Ly, Lz, ncells


# ========================================================================
def save_curve(plotnum, fname):
    """Save curve data"""
    SetActivePlots(plotnum)
    HideActivePlots()
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.outputToCurrentDirectory = 1
    SaveWindowAtts.outputDirectory = "."
    SaveWindowAtts.fileName = fname
    SaveWindowAtts.family = 0
    SaveWindowAtts.format = SaveWindowAtts.CURVE
    SaveWindowAtts.width = 1024
    SaveWindowAtts.height = 1024
    SaveWindowAtts.screenCapture = 0
    SaveWindowAtts.saveTiled = 0
    SaveWindowAtts.quality = 80
    SaveWindowAtts.progressive = 0
    SaveWindowAtts.binary = 0
    SaveWindowAtts.stereo = 0
    SaveWindowAtts.compression = SaveWindowAtts.None
    SaveWindowAtts.forceMerge = 0
    # NoConstraint, EqualWidthHeight, ScreenProportions
    SaveWindowAtts.resConstraint = SaveWindowAtts.ScreenProportions
    SaveWindowAtts.advancedMultiWindowSave = 0
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()
    HideActivePlots()


# ========================================================================
def get_lineout(field, datadir, p0, p1, steps=None, pfx=''):
    ret = AddPlot("Curve", "operators/Lineout/" + field, 1, 1)
    if ret == 0:
        return 0
    LineoutAtts = LineoutAttributes()
    LineoutAtts.point1 = tuple(truncate(i, 10) for i in p0)
    LineoutAtts.point2 = tuple(truncate(i, 10) for i in p1)
    LineoutAtts.interactive = 0
    LineoutAtts.ignoreGlobal = 0
    LineoutAtts.samplingOn = 0
    LineoutAtts.reflineLabels = 0
    SetOperatorOptions(LineoutAtts, 1)
    DrawPlots()
    SetActivePlots(0)
    HideActivePlots()

    # Loop over time
    m = GetMetaData(GetWindowInformation().activeSource)
    if steps is None:
        steps = range(TimeSliderGetNStates())
    for tidx in steps:
        TimeSliderSetState(tidx)
        print(
            "Saving lineout at time : {0:f} (cycle: {1:d})".format(
                m.times[tidx], m.cycles[tidx]
            )
        )
        save_curve(
            0,
            os.path.join(
                ".", datadir, pfx + field + "_{0:d}_{1:.16f}".format(tidx, m.times[tidx])
            ),
        )

    DeleteAllPlots()

    return 1


# ========================================================================
#
# Main
#
# ========================================================================

# Timer
start = time.time()

# Parse arguments
parser = argparse.ArgumentParser(description="A simple post-processing tool")
parser.add_argument(
    "-f", "--folder", help="Folder where plot files are stored", type=str, default="."
)
parser.add_argument("-i", "--iname", help="Input file name", type=str, required=True)
args, unknown = parser.parse_known_args()

# Get the list of data
if not os.path.exists(args.folder):
    print("Invalid folder as input argument (does not exist)")
    sys.exit()
return_code = sp.call(
    "ls -1v {0:s}/plt*/Header | tee movie.visit".format(args.folder), shell=True
)

# Store in directory
datadir = os.path.join(args.folder, "data")
if not os.path.exists(datadir):
    os.makedirs(datadir)

# Simulation attributes
icname = "ic.txt"
Lx, Ly, Lz, ncells = get_simulation_params(icname, args.iname)
dx = Lx / ncells[0]

# Open files
OpenDatabase("localhost:movie.visit", 0)

# Define expressions
DefineScalarExpression("magvel2", "sqr(magvel)")
DefineScalarExpression("cs_mod", "sqrt(pressure/density)")
DefineScalarExpression(
    "sij2",
    "sqr(gradient(x_velocity)[0]) + sqr(gradient(y_velocity)[1]) + sqr(gradient(z_velocity)[2]) + 2.0 * sqr(0.5 * (gradient(x_velocity)[1] + gradient(y_velocity)[0])) + 2.0 * sqr(0.5 * (gradient(x_velocity)[2] + gradient(z_velocity)[0])) + 2.0 * sqr(0.5 * (gradient(y_velocity)[2] + gradient(z_velocity)[1]))",
)
DefineScalarExpression("u2", "sqr(x_velocity)")

# Integrate on the whole domain
AddPlot("Pseudocolor", "magvel2", 1, 1)
AddPlot("Pseudocolor", "cs_mod", 1, 1)
AddPlot("Pseudocolor", "sij2", 1, 1)
AddPlot("Pseudocolor", "u2", 1, 1)
AddPlot("Pseudocolor", "Temp", 1, 1)
DrawPlots()
SetQueryFloatFormat("%g")
QueryOverTimeAtts = GetQueryOverTimeAttributes()
QueryOverTimeAtts.timeType = QueryOverTimeAtts.DTime  # Cycle, DTime, Timestep
QueryOverTimeAtts.startTimeFlag = 0
QueryOverTimeAtts.startTime = 0
QueryOverTimeAtts.endTimeFlag = 0
QueryOverTimeAtts.endTime = 1
QueryOverTimeAtts.strideFlag = 0
QueryOverTimeAtts.stride = 1
QueryOverTimeAtts.createWindow = 1
QueryOverTimeAtts.windowId = 2
SetQueryOverTimeAttributes(QueryOverTimeAtts)
SetActivePlots(0)
QueryOverTime("Weighted Variable Sum", do_time=1)
SetActivePlots(1)
QueryOverTime("Weighted Variable Sum", do_time=1)
SetActivePlots(2)
QueryOverTime("Weighted Variable Sum", do_time=1)
SetActivePlots(3)
QueryOverTime("Weighted Variable Sum", do_time=1)
SetActivePlots(4)
QueryOverTime("Max", do_time=1)

# Hide all curve plots
SetActiveWindow(2)
SetActivePlots(0)
HideActivePlots()
SetActivePlots(1)
HideActivePlots()
SetActivePlots(2)
HideActivePlots()
SetActivePlots(3)
HideActivePlots()
SetActivePlots(4)
HideActivePlots()

# Save
save_curve(0, os.path.join(".", datadir, "magvel2"))
save_curve(1, os.path.join(".", datadir, "cs_mod"))
save_curve(2, os.path.join(".", datadir, "Sij2"))
save_curve(3, os.path.join(".", datadir, "u2"))
save_curve(4, os.path.join(".", datadir, "maxtemp"))

# Clear everything
DeleteAllPlots()

# Lineouts through time
field = "x_velocity"
p0 = (0, 0.5 * Ly, 0.5 * Lz)
p1 = (Lx + 1e-3, 0.5 * Ly, 0.5 * Lz)
get_lineout(field, datadir, p0, p1)

field = "y_velocity"
p0 = (0.5 * Lx, 0, 0.5 * Lz)
p1 = (0.5 * Lx, Ly, 0.5 * Lz)
get_lineout(field, datadir, p0, p1)

field = "z_vorticity"
p0 = (0, 0.5 * Ly, 0.5 * Lz)
p1 = (Lx, 0.5 * Ly, 0.5 * Lz)
get_lineout(field, datadir, p0, p1)

field = "Temp"
p0 = (0.5 * Lx, 0, 0.5 * Lz)
p1 = (0.5 * Lx, Ly, 0.5 * Lz)
get_lineout(field, datadir, p0, p1)

field = "Y_lp_H2_rp_"
p0 = (0.5 * Lx, 0, 0.5 * Lz)
p1 = (0.5 * Lx, Ly, 0.5 * Lz)
get_lineout(field, datadir, p0, p1)

field = "Y_lp_O2_rp_"
p0 = (0.5 * Lx, 0, 0.5 * Lz)
p1 = (0.5 * Lx, Ly, 0.5 * Lz)
get_lineout(field, datadir, p0, p1)

field = "Y_lp_OH_rp_"
p0 = (0.5 * Lx, 0, 0.5 * Lz)
p1 = (0.5 * Lx, Ly, 0.5 * Lz)
get_lineout(field, datadir, p0, p1)

field = "heat_release"
p0 = (0.5 * Lx, 0, 0.5 * Lz)
p1 = (0.5 * Lx, Ly, 0.5 * Lz)
get_lineout(field, datadir, p0, p1)

# initial lineouts
field = "Temp"
p0 = (0, 0.5 * Ly, 0.5 * Lz)
p1 = (Lx, 0.5 * Ly, 0.5 * Lz)
get_lineout(field, datadir, p0, p1, steps=[0], pfx='init_')

field = "Y_lp_H2_rp_"
p0 = (0, 0.5 * Ly, 0.5 * Lz)
p1 = (Lx, 0.5 * Ly, 0.5 * Lz)
get_lineout(field, datadir, p0, p1, steps=[0], pfx='init_')

field = "Y_lp_O2_rp_"
p0 = (0, 0.5 * Ly, 0.5 * Lz)
p1 = (Lx, 0.5 * Ly, 0.5 * Lz)
get_lineout(field, datadir, p0, p1, steps=[0], pfx='init_')

# output timer
end = time.time() - start
print("Elapsed time " + str(timedelta(seconds=end)) + " (or {0:f} seconds)".format(end))

sys.exit()
