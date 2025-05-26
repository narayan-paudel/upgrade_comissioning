#!/usr/bin/env python

import os
import glob

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

from pathlib import Path
home = str(Path.home())
plotFolder = home+"/research_ua/icecube/upgrade/upgrade_comissioning/plots"

icetray_src = "/Users/epaudel/research_ua/icecube/software/icetray/src/"
tilt_file = icetray_src + "ice-models/resources/models/PPCTABLES/misc/cx.dat.0"


with open(tilt_file,"r") as fh:
    lines = fh.readlines()
    strings = []
    oms = []
    nxs = []
    nys = []
    nzs = []
    errors = []
    for x in lines:
        strings.append(x.split(" ")[0])
        oms.append(x.split(" ")[1])
        nxs.append(float(x.split(" ")[2]))
        nys.append(float(x.split(" ")[3]))
        nzs.append(float(x.split(" ")[4]))
        errors.append(float(x.split(" ")[5]))
# print(sorted(list(set(strings))))
# print(sorted(list(set(oms))))

nz_unit = np.array([0,0,-1])
first_tilt_unit = np.array([nxs[0],nys[0],nzs[0]])/np.sqrt(nxs[0]**2+nys[0]**2+nzs[0]**2)
print(nz_unit,first_tilt_unit)
print(np.linalg.norm(first_tilt_unit))

def to_spherical(x,y,z):
    '''
    converts cartesian coordinates (x,y,z) to spherical (r,theta,phi)
    '''
    r = np.sqrt(x*x+y*y+z*z)
    theta = np.arctan2(np.sqrt(x*x + y*y),z)
    phi = np.arctan2(y,x)
    return r,theta,phi

def angle_between_vectors(vector1,vector2):
    '''
    calculates angle between two vectors
    '''
    v1 = np.asarray(vector1)
    v2 = np.asarray(vector2)
    return np.arccos(np.clip(np.dot(v1,v2),-1.0,1.0))

tilt_angles = []
for inx,iny,inz in zip(nxs,nys,nzs):
    tilt_angles.append(np.arctan2(np.sqrt(inx*inx + iny*iny),inz))

print(angle_between_vectors(nz_unit,first_tilt_unit)*180/np.pi)
r,theta,phi = to_spherical(first_tilt_unit[0],first_tilt_unit[1],first_tilt_unit[2])
print(np.rad2deg(theta),np.rad2deg(phi))

print(nxs[0],nys[0],nzs[0])
print(first_tilt_unit)

def tilt_hist(tilt_angles):
    tilt_angles = [np.rad2deg(itilt) for itilt in tilt_angles]
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    bins = np.linspace(-1,60,122)
    ax.hist(tilt_angles,bins=bins,histtype="step",color="#a6cee3",linewidth=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r"tilt angle $\theta^{\circ}$", fontsize=22)
    ax.set_ylabel("count", fontsize=22)
    # ax.set_xlim(0,100)
    # ax.set_ylim(0.9,5*10**3)
    ax.set_yscale("log")
    ax.grid(True,alpha=0.6)
    # ax.legend(fontsize=14)
    # ax.legend(fontsize=8,ncols=2,bbox_to_anchor=(0.5, 1.05),loc="center")
    plt.savefig(plotFolder+f"/tilt_dist.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/tilt_dist.pdf",transparent=False,bbox_inches='tight')
    plt.close()
tilt_hist(tilt_angles)
print(max(tilt_angles)*180/np.pi)
tilt_angles = [np.rad2deg(itilt) for itilt in tilt_angles]
for itilt,istring,idom in zip(tilt_angles,strings,oms):
    if itilt > 30:
        print(f"String {istring} DOM {idom} tilt {itilt:.1f} degrees")