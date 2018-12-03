#
# Post-process data using VisIt. Generates data files of
# - the spatial average of u^+v^2+w^2
# - the spatial average of sound speed
# - the spatial average of u^2
# - the spatial average of dudx^2
#
# Usage:
#    visit -nowin -cli -s /path/to/pp_aux_vars.py
#

# ========================================================================
#
# Imports
#
# ========================================================================
import sys
import subprocess as sp


# ========================================================================
#
# Function definitions
#
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
#
# Main
#
# ========================================================================

# Get the list of data
return_code = sp.call('ls -1 plt*/Header | tee movie.visit', shell=True)

# Open files
OpenDatabase("localhost:movie.visit", 0)

# Define expressions
DefineScalarExpression("magvel2", "sqr(magvel)")
DefineScalarExpression("cs_mod", "sqrt(pressure/density)")
DefineScalarExpression("dudx2", "sqr(gradient(x_velocity)[0])")
DefineScalarExpression("u2", "sqr(x_velocity)")

# Integrate on the whole domain
AddPlot("Pseudocolor", "magvel2", 1, 1)
AddPlot("Pseudocolor", "cs_mod", 1, 1)
AddPlot("Pseudocolor", "dudx2", 1, 1)
AddPlot("Pseudocolor", "u2", 1, 1)
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

# Save the square of the magnitude of velocity
save_curve(0, "magvel2")

# Save the modified sound speed data
save_curve(1, "cs_mod")

# Save the square of the x_velocity gradient
save_curve(2, "dudx2")

# Save the square of the x_velocity
save_curve(3, "u2")

sys.exit()
