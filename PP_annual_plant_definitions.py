"""level plots of the different definitions for fitness and niche
differences for the annual plant model"""

import matplotlib.pyplot as plt
from matplotlib.cm import rainbow
import matplotlib.patches as patches
import numpy as np

import plot_Figure1 as apd

def save_ND(counter = [1]):
    print(counter[0])
    plt.savefig("PP_slides/ND in annual plants_{}.png".format(counter[0]),
                transparent = "True")
    counter[0] += 1
    

###############################################################################
# plotting the results for ND
#plt.style.use('dark_background')
fig = plt.figure(figsize = (8,8))

fs_label = 20
fs_axis = 18
# axis labels
plt.xlabel(r'Interspecific interaction ($\alpha$)', fontsize = fs_label)
plt.ylabel(r'Niche difference $(\mathcal{N})$', fontsize = fs_label)
ND_range = [-0.5,1.5]
plt.axis([min(apd.interspec), max(apd.interspec)] + ND_range)

plt.xticks([])
plt.yticks([])

save_ND()

# interspecific competition absent
plt.plot(0,1, 'o', color = "k", markersize = 10)
plt.xticks([0], fontsize = fs_axis)
plt.yticks([1], fontsize = fs_axis)
save_ND()

# interspecific equal to intraspecific
plt.plot(apd.sign*1,0, 'o', color = "k", markersize = 10)
plt.xticks([apd.sign*1, 0], fontsize = fs_axis)
plt.yticks([1, 0], fontsize = fs_axis)
save_ND()

offset_y = 0.15
ax = plt.gca()

# normal competition
rect_norm = patches.Rectangle([0,0],apd.sign,1, fill = False)
ax.add_patch(rect_norm)
plt.text(apd.sign*1/2,ND_range[0]+offset_y,
         "negative,\nweaker than\nintraspecific", 
         ha = "center",va = "center", fontsize = fs_axis)
save_ND()

# positive interactions
rect_facilitation = patches.Rectangle([apd.interspec[0],1],
            -apd.interspec[0], ND_range[1]-1, fill = False, linestyle = ":")
ax.add_patch(rect_facilitation)
plt.text(apd.interspec[0]/2,ND_range[0]+offset_y,"positive", 
         ha = "center", fontsize = fs_axis)
save_ND()

# stronger than intraspecific competition
plt.text((apd.sign*1+apd.interspec[-1])/2-0.09,ND_range[0]+offset_y,
          "negative,\nstronger than\nintraspecific",
         ha = "center",va = "center", fontsize = fs_axis)
rect_comp = patches.Rectangle([apd.sign*1,0],apd.interspec[-1], ND_range[0]
                              , fill = False, linestyle = "--")
ax.add_patch(rect_comp)
save_ND()

keys = ["Chesson (2003)", "Zhao et al. (2016)", "Bimler et al. (2018)",
        "Adler et al. (2007)", "Carmel et al. (2017)",
        "Carroll et al. (2011)", "Saavedra et al. (2017)", 
        "Godoy & Levine (2014)", "Spaak & De Laender"]
colors =  {keys[i]: rainbow(np.linspace(0, 1, len(keys)))[i]
                for i in range(len(keys))}

lw_new = 4
lw_old = 2
# plot definitions that don't hit any point at all               
for key in keys[:5]:
    plt.plot(apd.interspec, apd.ND[key], label = key, linewidth = lw_new, 
             color = colors[key])
plt.legend(loc = "upper left", framealpha = 0)
save_ND()
# reduce contrast of already shown definitions                
for l in ax.lines:
    l.set_alpha(0.3)
    l.set_linewidth(lw_old)
    
# plot definitions that hit one point                
for key in keys[5:7]:
    plt.plot(apd.interspec, apd.ND[key], label = key, linewidth = lw_new, 
             color = colors[key])
plt.legend(loc = "upper left", framealpha = 0)
save_ND()
# reduce contrast of already shown definitions                
for l in ax.lines:
    l.set_alpha(0.3)
    l.set_linewidth(lw_old)
    
# plot definitions that hit both point
key = keys[7]            
plt.plot(apd.interspec, apd.ND[key], label = key, linewidth = lw_new,
             color = colors[key])
plt.legend(loc = "upper left", framealpha = 0)
save_ND()
# reduce contrast of already shown definitions                
for l in ax.lines:
    l.set_alpha(0.3)
    l.set_linewidth(lw_old)

key = keys[-1]              
plt.plot(apd.interspec, apd.ND[key], label = key, linewidth = lw_new,
             color = colors[key])
plt.legend(loc = "upper left", framealpha = 0)
save_ND()