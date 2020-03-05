"plot the emission spectrum of the incoming light"

import matplotlib.pyplot as plt
import pandas as pd

k_spec = pd.read_csv("Lights_and_absorptions.csv")
a = pd.read_csv("26.csv", skiprows = 1, skipfooter = 10)
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
ax_I_in.bar(k_spec["lambda"][::dist]+2*dist, k_spec["I_in"][::dist],
            width = dist,
           color = "red")
plt.ylabel(r"Incoming light [mumol photons $m^{-2}s^{-1}$]")

fig.savefig("AP_figure_incoming_light.pdf")