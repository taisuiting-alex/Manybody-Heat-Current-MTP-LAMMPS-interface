from pylab import *

aw = 1.5
fs = 24
lw = 2
font = {'size'   : fs,'family':'Times New Roman'}
matplotlib.rc('font', **font)
matplotlib.rc('axes' , linewidth=aw)


#Generate setting:
cut=500
dt = 0.001  # ps 
Ns = 1000  # Sample interval
BLOCK_VOLUME=24.6*3.35*4260 #A**3
BLOCK_LENGTH = 4260 #A
CROSS_AREA=BLOCK_VOLUME/BLOCK_LENGTH #A**2

def lammps_heatflux(path):
    thermo = loadtxt(path + "/compute_Energy_Temp.out")
    jp = loadtxt(path + "/compute_HeatFlux.out")

    Ein = thermo[cut:, 1]
    Eout = thermo[cut:, 2]
    Etol = (Eout - Ein) / 2 / 1000 #in units of KeV
    Etol = Etol - Etol[0]
    
    #t = dt * np.arange(1, len(Etol) + 1) * Ns / 1000 #unit in ns #fix scale
    tt= (thermo[cut:,0]-thermo[cut,0])/1000*dt #unit in ns # adjust by dt
    
    jpy = jp[cut:, 2] - jp[cut:, 5]
    jpy = jpy*-1            ##negative y direction
    jpy = jpy / BLOCK_VOLUME *CROSS_AREA * 1000 #in units of eV/ns=eVA/ps*(1/A**3)*A**2*1000ps/ns
    accum_jpy = cumsum(jpy) * dt / 1000 #in units of KeV=eV/ns*dt*ns/1000ps*KeV/1000eV
    
    return tt, accum_jpy, Etol


import sys
import numpy as np
import matplotlib.pyplot as plt
directory_list=sys.argv[1:]
n_plot=len(directory_list)
plt.figure(1,figsize=(12, 6*int(n_plot/2+n_plot%2)))
plt.figure(2,figsize=(12, 6*int(n_plot/2+n_plot%2)))
for i in range(n_plot):
    print("Reading : "+directory_list[i])
    plt.figure(1)
    plt.subplot(int(n_plot/2+n_plot%2), 2, i+1)
    
    tl = 6
    tw = 1.5
    tlm = 3    
    plt.gca().tick_params(which='major', length=tl, width=tw)
    plt.gca().tick_params(which='minor', length=tlm, width=tw)
    plt.gca().tick_params(axis='both', direction='out', right=False, top=False)
    
    t_mtp, jp_mtp, etol_mtp = lammps_heatflux(directory_list[i])
    
    plt.plot(t_mtp,jp_mtp, lw = lw, ls = "-", label = r"MTP  Heat Current, "+'$J_{pot}$')
    plt.plot(t_mtp, etol_mtp, lw = lw, ls = "--", label = r"MD Langervin Thermostat Heat Current, "+'$\\frac{(E_{in} + E_{out})}{2}$')
    plt.legend(fontsize=18,loc=2, bbox_to_anchor=(-0.30, -0.7, 0.5,0.5),frameon=False,framealpha=1.0)
    
    left=2.0-cut/1000
    plt.xlim([0, left])
    plt.gca().set_xticks(linspace(0, left, 4))
    
    m=1.8
    ylim([0, m])
    plt.gca().set_yticks(linspace(0, m, 4))
    plt.xlabel('Time (ns)')#,fontsize=16)
    plt.ylabel('Accumulated heat (keV)')#,fontsize=16)
    #plt.title("(?) "+directory_list[i])
    


plt.figure(1)
subplots_adjust(wspace = 0.3, hspace = 0.35)
savefig("figure_out.png", bbox_inches='tight')


