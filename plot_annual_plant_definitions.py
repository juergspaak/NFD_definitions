"""level plots of the different definitions for fitness and niche
differences for the annual plant model"""

import matplotlib.pyplot as plt
from matplotlib.cm import rainbow
import matplotlib.patches as patches
import numpy as np
from scipy.optimize import brentq
plt.rcParams["font.family"] = 'Times New Roman'
np.seterr(all = "ignore")

rep = 5000

A11 = 1
A22 = 1
A21 = 0.7*np.ones(rep)

lamb1 = 1.5
lamb2 = 3

A12 = np.linspace(-A22/(lamb2-1)+1e-3,2.5, rep)
sign = -1
interspec = sign*A12*A21/A11/A22 # interspecific interaction strength

ND = {}
FD = {}
ND_bound = {}
FD_bound = {}

# Definition according to Adler et al. 2007
key = "Adler et al. (2007)"
ND[key] = np.log(lamb2/(1+A12/A22*(lamb2-1)))
FD[key] = np.log(lamb1/lamb2*np.ones(A12.shape))

# Definition according to Godoy et al.
key = "Godoy & Levine (2014)"
ND[key] = 1-np.sqrt(A12*A21/(A22*A11))
FD[key] = (lamb1-1)/(lamb2-1)*np.sqrt(A21*A22/(A11*A12))

# Definitions by experimental methods
# zero growth rate
f_i_0 = np.log([lamb1,lamb2])

# invasion growth rate
r_i = np.array([np.log(lamb1/(1+A12/A22*(lamb2-1))),
                np.log(lamb2/(1+A21/A11*(lamb1-1)))])

# Definition according to Carroll et al.
S_i = (f_i_0.reshape(-1,1)-r_i)/f_i_0.reshape(-1,1)
key = "Carroll et al. (2011)"
ND[key] = 1-np.exp(np.mean(np.log(S_i),axis = 0))
FD[key] = np.exp(np.std(np.log(S_i),axis = 0))

# Definition according to Zhao et al.
key = "Zhao et al. (2016)"
ND[key] = 1+np.sum(r_i, axis = 0)
FD[key] = (np.log((lamb1-1)/A11)-np.log((lamb2-1)/A22))*np.ones(A12.shape)
    
# Definition according to Chesson
key = "Chesson (2003)"
phi_i = np.array([A11/lamb1,A22/lamb2])
ND[key] = 0.5*np.sum(r_i/phi_i.reshape(-1,1),axis = 0)
FD[key] = (r_i/phi_i.reshape(-1,1)-ND[key])[0]

# Definition according to Bimler
key = "Bimler et al. (2018)"

# change to alpha' according to Bimler et al
A11_ = A11/(lamb1-1)
A12_ = A12/(lamb1-1)
A21_ = A21/(lamb2-1)
A22_ = A22/(lamb2-1)

ND[key] = 1 - np.exp((A12_-A11_)+(A21_-A22_))
FD[key] = np.exp(-(A12_+A11_)+(A21_+A22_))

# Definition according to Saavedra
key = "Saavedra et al. (2017)"

ND[key] = 2/np.pi*np.arcsin((A11*A22-A12*A21)/
                    (np.sqrt((A11**2+A21**2)*(A12**2+A22**2))))
r = np.array([[lamb1-1], [lamb2-1]])
rc = 0.5*(1/np.sqrt(A11**2+A21**2)*np.array([A11*np.ones(rep),A21])+
          1/np.sqrt(A12**2+A11**2)*np.array([A12,A22*np.ones(rep)]))
FD[key] = 180/np.pi*np.arccos(np.sum(r*rc,axis = 0)/np.linalg.norm(r)
                                        /np.linalg.norm(rc,axis = 0))

# Definition according to Karmel
key = "Carmel et al. (2017)"
# mortality rate, assumed to be minimal invasion growth rate of all
mort = -2*np.nanmin(r_i)
geo_mean = np.sqrt(np.prod(r_i/mort+1, axis = 0))
NO_1 = -(1-2*geo_mean**2) + np.sqrt((1-2*geo_mean**2)**2-1)
NO_2 = -(1-2*geo_mean**2) - np.sqrt((1-2*geo_mean**2)**2-1)


ND[key] = 1 - np.nanmin([NO_1,NO_2], axis = 0)
# add minimal value before nans:
ND[key][np.argmax(np.isnan(ND[key]))] = 0
FD[key] = np.exp(np.std(np.log(r_i/mort+1),axis = 0))**2
s = r_i/mort +1
test1 = np.amin([s[0]/s[1], s[1]/s[0]], axis = 0)

# Definition accoding to Spaak


c = np.empty(rep)
denom = np.array([np.log(1+A12/A22*(lamb2-1)),
                  np.log(1+A21/A11*(lamb1-1))])
fac = np.array([A11/A22*(lamb2-1),A22/A11*(lamb1-1)])

def NO_spaak(c,denom):
    return np.array([denom[0]/np.log(1+c*fac[0]),denom[1]/np.log(1+fac[1]/c)])
def NO_equate(c,denom):
    NO = np.abs(NO_spaak(c,denom))
    return NO[0]-NO[1]

for l in range(rep):
    try:
        c[l] = brentq(NO_equate,0,1e3, args = (denom[:,l]))
    except ValueError:
        c[l] = brentq(NO_equate,0,1e10, args = (denom[:,l]))
key = "Spaak & De Laender"   
ND[key] = 1-NO_spaak(c, denom)[0]
FD[key] = np.log(lamb1/(1+c*A11/A22*(lamb2-1)))/np.log(lamb1)

keys = ["Chesson (2003)","Carroll et al. (2011)", "Zhao et al. (2016)",
        "Godoy & Levine (2014)", "Adler et al. (2007)", "Bimler et al. (2018)",
        "Carmel et al. (2017)", "Saavedra et al. (2017)", "Spaak & De Laender"]

colors =  {keys[i]: rainbow(np.linspace(0, 1, len(keys)))[i]
                for i in range(len(keys))}


if __name__ == "__main__":
    ###########################################################################
    # plotting the results for ND
    fig = plt.figure(figsize = (10,10))
    
    ND_range = [-0.5,1.5]
    
    rect_facilitation = patches.Rectangle([interspec[0],1],
                -interspec[0], ND_range[1]-1, fill = False, linestyle = ":")
    rect_norm = patches.Rectangle([0,0],sign,1, fill = False)
    rect_comp = patches.Rectangle([sign*1,0],interspec[-1], ND_range[0]
                                  , fill = False, linestyle = "--")
    ax = plt.gca()
    ax.add_patch(rect_norm)
    ax.add_patch(rect_facilitation)
    ax.add_patch(rect_comp)
              
    # plot NFD parameters          
    for key in keys:
        plt.plot(interspec, ND[key], label = key, linewidth = 2, alpha = 1, 
                 color = colors[key])
    
    # add black dots
    plt.plot(0,1, 'o', color = "black", markersize = 10)
    plt.plot(sign*1,0, '^', color = "black", markersize = 10)
    
    
    # layout
    plt.legend(loc = "upper left", fontsize = 12)
    
    # axis limits
    plt.xlim(min(interspec), max(interspec))
    plt.xticks([0,sign*1])
    plt.ylim(*ND_range)
    plt.yticks([0,1])
    
    # labeling of interaction
    fs = 16
    offset_y = 0.15
    plt.text(interspec[0]/2,ND_range[0]+offset_y,"positive", 
             ha = "center", fontsize = fs, backgroundcolor = "white")
    plt.text(sign*1/2,ND_range[0]+offset_y,
             "negative,\nweaker than\nintraspecific", 
             ha = "center",va = "center", fontsize = fs)
    plt.text((sign*1+interspec[-1])/2-0.09,ND_range[0]+offset_y,
              "negative,\nstronger than\nintraspecific",
             ha = "center",va = "center", fontsize = fs)
    
    
    
    # axis labels
    plt.xlabel(r'Interspecific competition ($\alpha$)', fontsize = 16)
    plt.ylabel(r'Niche difference $(\mathcal{N})$', fontsize = 16)
    
    fig.savefig("ND_in_annual_plants.eps")
    
    ###########################################################################
    # plotting the results for ND
    fig = plt.figure(figsize = (9,9))
    
    ND_range = [-0.5,1.5]
              
    # plot NFD parameters          
    for key in keys:
        plt.plot(interspec, FD[key], label = key, linewidth = 2, alpha = 1, 
                 color = colors[key])
    
    # layout
    plt.legend(loc = "upper left", fontsize = 12)
    
    # axis limits
    plt.xlim(min(interspec), max(interspec))
    plt.xticks([0,sign*1])
    plt.ylim(-3,3)
    
    
    plt.axhline(y=0, color = "black", linestyle = ":")
    plt.axvline(x=0, color = "black", linestyle = ":")
    
    # axis labels
    plt.xlabel(r'Interspecific competition ($\alpha$)', fontsize = 16)
    plt.ylabel(r'Fitness difference $(\mathcal{F})$', fontsize = 16)
    
    fig.savefig("FD_in_annual_plants.pdf")