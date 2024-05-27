#!/usr/bin/env python

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams["font.family"] = "Serif"
matplotlib.rcParams["font.size"] = 10
matplotlib.rcParams["axes.labelsize"] = 10
matplotlib.rcParams["xtick.labelsize"] = 10
matplotlib.rcParams["ytick.labelsize"] = 10
matplotlib.rcParams["legend.fontsize"] = 10

fig = plt.figure(facecolor="white")
ax = fig.add_subplot(1, 1, 1)
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel("qv")
ax.set_ylabel("H")
ax.set_title("pipeline")

x = np.array([0.9539999999999999E+01,0.9260000000000000E+01,0.8900000000000000E+01,0.8539999999999999E+01,0.8170000000000000E+01,0.7800000000000000E+01,0.7410000000000000E+01,0.7020000000000000E+01,0.6640000000000000E+01,0.6240000000000000E+01,0.5840000000000000E+01,0.5440000000000000E+01,0.5030000000000000E+01,0.4630000000000000E+01,0.4220000000000000E+01,0.3810000000000000E+01,0.2760000000000000E+01,0.1590000000000000E+01])
y = np.array([0.5995174407958984E+01,0.5664133548736572E+01,0.5249705314636230E+01,0.4831446170806885E+01,0.4309559345245361E+01,0.3885838270187378E+01,0.3462457418441772E+01,0.2932367086410522E+01,0.2804669380187988E+01,0.2676398754119873E+01,0.2543397903442383E+01,0.2405666589736938E+01,0.2366906881332397E+01,0.2117384195327759E+01,0.2068806648254395E+01,0.1913046240806580E+01,0.1708378672599792E+01,0.1476971149444580E+01])

ax.plot(x,y,"b-o",linewidth=2,markersize=5,label="need")

x = np.array([0.9539999999999999E+01,0.9260000000000000E+01,0.8900000000000000E+01,0.8539999999999999E+01,0.8170000000000000E+01,0.7800000000000000E+01,0.7410000000000000E+01,0.7020000000000000E+01,0.6640000000000000E+01,0.6240000000000000E+01,0.5840000000000000E+01,0.5440000000000000E+01,0.5030000000000000E+01,0.4630000000000000E+01,0.4220000000000000E+01,0.3810000000000000E+01,0.2760000000000000E+01,0.1590000000000000E+01])
y = np.array([0.6000808715820312E+01,0.6947523593902588E+01,0.8314241409301758E+01,0.1036411857604980E+02,0.1129964828491211E+02,0.1223994159698486E+02,0.1335399913787842E+02,0.1417295932769775E+02,0.1485116863250732E+02,0.1584093475341797E+02,0.1653754234313965E+02,0.1730471038818359E+02,0.1799187088012695E+02,0.1854858016967773E+02,0.1909412765502930E+02,0.1953249168395996E+02,0.2006701850891113E+02,0.2059832191467285E+02])

ax.plot(x,y,"r-o",linewidth=2,markersize=5,label="provide")

ax.legend(loc="best")

plt.savefig("pipeline.png")

