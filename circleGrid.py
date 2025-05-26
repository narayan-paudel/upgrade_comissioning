#!/usr/bin/env python

import os
import glob

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pathlib import Path
home = str(Path.home())

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = home+"/research_ua/icecube/upgrade/upgrade_comissioning/plots/"
# angles = np.linspace(0,360,36001)
angles = np.linspace(0,90,9001)
angles_rad = [np.deg2rad(itheta) for itheta in angles]
x = [np.cos(itheta) for itheta in angles_rad]
y = [np.sin(itheta) for itheta in angles_rad]


zero_line_x = np.linspace(0,1,100)
zero_line_y = [0 for ix in zero_line_x]

# angles_coarse = np.linspace(0,90,91)
angles_coarse = np.linspace(0,360,37)
print(angles_coarse)

def rotate(x,y,theta):
    theta_r = np.deg2rad(theta)
    return np.cos(theta_r), np.sin(theta_r)


    



def plot_circleGrid(x,y):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    ax.plot(x,y,"-",label="azimuth")
    for iangle in angles_coarse:
        x_rot, y_rot = rotate(1,0,iangle)
        # ax.plot([x_rot,-x_rot],[y_rot,-y_rot],"-",lw=1.0)
        ax.plot([x_rot,0],[y_rot,0],"-",lw=1.0)
    # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.set_xlabel(r" x", fontsize=22)
    ax.set_ylabel(r" y", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    ax.set_aspect('equal')
    # ax.legend()
    # ax.set_xlim(0,1)
    # ax.set_ylim(0,1)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.axis('off')
    plt.savefig(plotFolder+f"/circleGrid.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/circleGrid.pdf",transparent=False,bbox_inches='tight')
    plt.close()


plot_circleGrid(x,y)