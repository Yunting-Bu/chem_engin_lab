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
ax.set_xlabel("Re")
ax.set_ylabel("Lambda")
ax.set_title("Plot of Lambda~Re")

x = np.array([0.7577062377929688E+03,0.1136559326171875E+04,0.1515412475585938E+04,0.1894265502929688E+04,0.2273118652343750E+04,0.2651971435546875E+04,0.3030824951171875E+04,0.3409677734375000E+04,0.3788531005859375E+04,0.7577062011718750E+04,0.1136559277343750E+05,0.1515412402343750E+05,0.1894265429687500E+05,0.2273118554687500E+05,0.2651971679687500E+05])
y = np.array([0.1115186572074890E+01,0.8619800209999084E+00,0.7697210907936096E+00,0.7699636220932007E+00,0.9062433242797852E+00,0.7061636447906494E+00,0.5715511441230774E+00,0.5004171729087830E+00,0.4448830783367157E+00,0.4448830783367157E+00,0.3921562135219574E+00,0.3577601313591003E+00,0.3333657383918762E+00,0.3111435472965240E+00,0.2955799102783203E+00])

ax.plot(x,y,"b-o",linewidth=2,markersize=5,label="Lambda~Re")
ax.set_xscale("log")
ax.set_yscale("log")

ax.legend(loc="best")

plt.savefig("Lambda~Re2.png")

