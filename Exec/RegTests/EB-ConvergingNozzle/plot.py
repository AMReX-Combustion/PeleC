import numpy as np
import pandas as pd
from scipy.optimize import root
import matplotlib.pyplot as plt

# Create Mach number area relation
def AoverAstarFromMachNumber(MachNumber, gamma=1.4):
    AoverAstarSquared = (1 / MachNumber**2) * (
        2 / (gamma + 1) * (1 + (gamma - 1) / 2 * MachNumber**2)
    ) ** ((gamma + 1) / (gamma - 1))
    return np.sqrt(AoverAstarSquared)


def MachNumberFromAoverAstar(AoverAstar, gamma=1.4, MachNumberGuess=0.2):
    ZeroFunc = lambda x: AoverAstarFromMachNumber(x, gamma) - np.array(AoverAstar)
    guess = np.zeros(np.array(AoverAstar).shape) + MachNumberGuess
    return abs(root(ZeroFunc, guess).x)


def P0OverPFromMachNumber(MachNumber, gamma=1.4):
    return (1 + 0.5 * (gamma - 1) * MachNumber**2) ** (gamma / (gamma - 1))


def AofX(X, inletdiam, outletdiam, nozzlestart, nozzleend):
    diam = outletdiam + (inletdiam - outletdiam) * (np.array(X) - nozzleend) / (
        nozzlestart - nozzleend
    )
    diam = np.clip(diam, min(inletdiam, outletdiam), max(inletdiam, outletdiam))
    return np.pi * (diam / 2) ** 2


# Nominal Conditions
M = 0.22
indiam = 10.0
outdiam = 10.0 / np.sqrt(2.0)
nozstart = 4.5
nozlen = 9

# load results
results = pd.read_csv("run_01/plt24000.slice", header=2, delim_whitespace=True)
results.columns = list(results.columns[1:]) + ["N/A"]
results["Area"] = AofX(results["x"], indiam, outdiam, nozstart, nozstart + nozlen)

AoverAstar0 = AoverAstarFromMachNumber(M)
Astar = results["Area"][0] / AoverAstar0
results["AoverAstar"] = results["Area"] / Astar
results["MachTheory"] = MachNumberFromAoverAstar(results["AoverAstar"])
results["P0OverPTheory"] = P0OverPFromMachNumber(results["MachTheory"])
P0 = results["pressure"][0] * P0OverPFromMachNumber(results["MachNumber"][0])
results["P0OverP"] = P0 / results["pressure"]

AoAs = AoverAstarFromMachNumber(M)
Mout = MachNumberFromAoverAstar(AoAs)

plt.figure("MachNumber", figsize=(5.5, 3))
plt.clf()
plt.xlabel("Axial Distance (cm)")
plt.ylabel("Mach Number")
plt.plot(results["x"], results["MachNumber"], label="PeleC")
plt.plot(results["x"], results["MachTheory"], label="Quasi-1D Flow")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("MachNumber.png", dpi=150)
plt.close()

plt.figure("Pressure", figsize=(5.5, 3))
plt.clf()
plt.xlabel("Axial Distance (cm)")
plt.ylabel("$p_0/p$")
plt.plot(results["x"], results["P0OverP"], label="PeleC")
plt.plot(results["x"], results["P0OverPTheory"], label="Quasi-1D Flow")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("Pressure.png", dpi=150)
plt.close()

# Compare different boundary conditions
data_00 = pd.read_csv("run_00/extremalog", delim_whitespace=True)
data_01 = pd.read_csv("run_01/extremalog", delim_whitespace=True)
data_10 = pd.read_csv("run_10/extremalog", delim_whitespace=True)

plt.figure("BCCompare", figsize=(5.5, 4.5))
plt.clf()
plt.plot(
    data_10["time"],
    data_10["pressure-min"] / 1e6,
    "k-",
    label="Inlet: Fixed P, Outlet: FOExtrap",
)
plt.plot(data_10["time"], data_10["pressure-max"] / 1e6, "k-")
plt.plot(
    data_00["time"],
    data_00["pressure-min"] / 1e6,
    "r-",
    label="Inlet: Interior P, Outlet: FOExtrap",
)
plt.plot(data_00["time"], data_00["pressure-max"] / 1e6, "r-")
plt.plot(
    data_01["time"],
    data_01["pressure-min"] / 1e6,
    "b-",
    label="Inlet: Interior P, Outlet: Characteristic",
)
plt.plot(data_01["time"], data_01["pressure-max"] / 1e6, "b-")
plt.xlabel("Time (s)")
plt.ylabel("Min/Max Pressure ($10^6$ dyn/cm$^2$)")
plt.legend(frameon=False)
plt.xlim([0, 0.024])
plt.ylim([0, 2.7])
plt.tight_layout()
plt.savefig("BCCompare.png", dpi=150)
plt.close()
