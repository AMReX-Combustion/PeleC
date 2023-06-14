import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df_deflt = pd.read_csv("data_vfracdef")
df_nosvf = pd.read_csv("data_vfrac0")
df_smcfl = pd.read_csv("data_vfracdef_smalldt")

plt.figure()

plt.plot(df_deflt["time"], df_deflt["min T"], "r-.")
plt.plot(df_deflt["time"], df_deflt["max T"], "r-.")
plt.plot(df_deflt["time"], df_deflt["avg T"], "r-", label="default")

plt.plot(df_smcfl["time"][::2], df_smcfl["min T"][::2], "kx")
plt.plot(df_smcfl["time"][::2], df_smcfl["max T"][::2], "kx")
plt.plot(df_smcfl["time"][::2], df_smcfl["avg T"][::2], "ko", label="small dt")

plt.plot(df_nosvf["time"], df_nosvf["min T"], "b-.")
plt.plot(df_nosvf["time"], df_nosvf["max T"], "b-.")
plt.plot(df_nosvf["time"], df_nosvf["avg T"], "b-", label="no small vfrac")

plt.legend()
plt.xlabel("Time")
plt.ylabel("Min/Max/Mean Temperature")

plt.savefig("comparison.png")

plt.figure()
plt.yscale("log")
plt.plot(df_deflt["time"], np.abs(400 - df_deflt["avg T"]) / 100, "r-", label="default")
plt.plot(
    df_smcfl["time"][::2],
    np.abs(400 - df_smcfl["avg T"][::2]) / 100,
    "ko",
    label="small dt",
)
plt.plot(
    df_nosvf["time"],
    np.abs(400 - df_nosvf["avg T"]) / 100,
    "b-",
    label="no small vfrac",
)

plt.legend()
plt.xlabel("Time")
plt.ylabel("abs(T-Twall)/T0-Twall")
plt.savefig("comparison-log.png")
