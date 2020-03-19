"""
@author: J.W.Spaak
Create Fiugre A3
"""

import matplotlib.pyplot as plt
import pandas as pd

k_spec = pd.read_csv("Lights_and_absorptions.csv")
n_bars = 15
dist = len(k_spec)//n_bars

# change units from absorption spectrum, current unit: mul/cells/cm
# change to ml/cells/m
unit_conv = 1000*100
fig = plt.figure()
plt.bar(k_spec["lambda"][::dist], unit_conv*(k_spec["BS4; 0"])[::dist],
           width = dist,
           color = "black")
plt.bar(k_spec["lambda"][::dist]+dist, unit_conv*(k_spec["BS5; 0"])[::dist],
           width = dist,
           color = "lightgrey")

plt.xlabel(r"wavelength $[nm]$")
plt.ylabel(r"absorption $[ml\cdot cells^{-1}\cdot m^{-1}]$")

plt.yticks([0,0.3,0.6])

ax_I_in = plt.gca().twinx()
ax_I_in.plot(k_spec["lambda"], k_spec["I_in"],
           color = "red")
plt.yticks([0,0.5,1])
plt.ylim([0,None])
plt.ylabel(r"Incoming light [mumol photons $m^{-2}s^{-1}$]")

fig.savefig("FigureA3.pdf")