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
# print(degg_list)

colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']


magnetometer_meas = np.asarray([[-9.924e-06, -4.78807e-05, -4.18883e-05],
[-9.51476e-06, -4.696e-05, -4.16837e-05],
[-9.19322e-06, -4.76761e-05, -4.16691e-05],
[-9.41245e-06, -4.79684e-05, -4.28822e-05],
[-9.70476e-06, -4.79538e-05, -4.16106e-05],
[-9.29553e-06, -4.81e-05, -4.27653e-05],
[-9.19322e-06, -4.82315e-05, -4.24876e-05],
[-9.60246e-06, -4.78807e-05, -4.22099e-05],
[-1.01286e-05, -4.753e-05, -4.26191e-05],
[-9.47091e-06, -4.80561e-05, -4.21953e-05]])#mDAB facing window, in T
magnetometer_meas_T = magnetometer_meas.T
mx = magnetometer_meas_T[0]
my = magnetometer_meas_T[1]
mz = magnetometer_meas_T[2]
mx_mean = np.mean(mx)
# print(np.mean(mx),np.std(mx))


meas_1 = np.asarray([[-1.49956e-05, -4.72084e-05, -4.57322e-05],
[2.91728e-05, -2.54165e-05, -4.0076e-05],
[3.50775e-06, 1.21602e-05, -4.25899e-05],
[-2.93774e-05, -9.76323e-06, -4.43145e-05],
[-1.49956e-05, -4.72084e-05, -4.57322e-05]])*10**6 #[x,y,z] in uT of mDAB2 facing window,PC, Chamber, Powersupply 

meas_2 = np.asarray([[-9.98246e-06, -4.77638e-05, -4.18006e-05],
[2.97866e-05, -2.34581e-05, -3.99445e-05],
[5.247e-06, 1.28764e-05, -4.24291e-05],
[-3.05759e-05, -9.25168e-06, -4.43876e-05],
[-9.98246e-06, -4.77638e-05, -4.18006e-05]])*10**6

meas_3 = np.asarray([[-6.7378e-06, -4.886e-05, -4.22245e-05], 
[2.94358e-05, -2.56504e-05, -3.98129e-05],
[4.22391e-06, 1.24671e-05, -4.13914e-05],
[-3.05759e-05, -9.55861e-06, -4.40368e-05]
,[-6.7378e-06, -4.886e-05, -4.22245e-05]])*10**6

def plot_magnet_field():
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    angles = [0,90,180,270,360]
    for imeas in [meas_1,meas_2,meas_3]:
        ax.plot(angles,imeas.T[0],"--o",color="#66c2a5",label="B$_{x}$",lw=2.5,alpha=1)
        ax.plot(angles,imeas.T[1],"--o",color="#fc8d62",label="B$_{y}$",lw=2.5,alpha=1)
        ax.plot(angles,imeas.T[2],"--o",color="#8da0cb",label="B$_{z}$",lw=2.5,alpha=1)
        ax.plot(angles,np.sqrt(imeas.T[0]**2+imeas.T[1]**2+imeas.T[2]**2),"--o",color="gray",label="B",lw=2.5,alpha=1)

    ax.set_xticks([0,90,180,270,360])
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
    ax.set_ylabel(r" B [$\mu$T]", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=1,bbox_to_anchor=(0.98, 0.70),loc="right")
    plt.savefig(plotFolder+f"/magnetometerxyz.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometerxyz.pdf",transparent=False,bbox_inches='tight')
    plt.close()
# plot_magnet_field()

def plot_magnet_field_polar():
    fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=1,ncols=1)
    # ax = fig.add_subplot(gs[0])
    ax = fig.add_subplot(polar=True)
    x_angles = [0,np.pi/2,np.pi,3*np.pi/2,2*np.pi]
    for imeas in [meas_1,meas_2,meas_3]:
        ax.plot(x_angles,imeas.T[0],"--o",color="#66c2a5",label="B$_{x}$",lw=2.5,alpha=1)
        ax.plot(x_angles,imeas.T[1],"--o",color="#fc8d62",label="B$_{y}$",lw=2.5,alpha=1)
        ax.plot(x_angles,imeas.T[2],"--o",color="#8da0cb",label="B$_{z}$",lw=2.5,alpha=1)
        ax.plot(x_angles,np.sqrt(imeas.T[0]**2+imeas.T[1]**2+imeas.T[2]**2),"--o",color="gray",label="B",lw=2.5,alpha=1)

    ax.set_rticks([-40,0,40])
    ax.set_rmax(70)
    # ax.set_rticks([0.5, 1, 1.5, 2])  # Less radial ticks
    # ax.set_rlabel_position(48)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    # ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=16)
    # ax.set_rlabel_position(45)
    ax.annotate(r" B [$\mu$T]",xy=[0,0.07], rotation=20,
            ha="left", va="bottom")
    # ax.set_rlabel(r" B [$\mu$T]", rotation=45, size=11)
    # ax.set_ylabel(r" B [$\mu$T]", fontsize=16)
    ax.grid(True,alpha=0.6)
    ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(1.05, 1.12),loc="right")
    plt.savefig(plotFolder+f"/magnetometerxyzPolar.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometerxyzPolar.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_magnet_field_polar()


magnetometer_readings = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/magnetometer_reading/magnetometer_measurement.xlsx"
import pandas as pd


df = pd.read_excel(magnetometer_readings,header=0)
rotations = df.columns.values
# print(df.columns.values)

def tilt_angle(nx,ny,nz):
    '''
    calculate tilt angle
    '''
    return np.arctan2(np.sqrt(nx*nx + ny*ny),nz)

def to_spherical(x,y,z):
    '''
    converts cartesian coordinates (x,y,z) to spherical (r,theta,phi)
    '''
    r = np.sqrt(x*x+y*y+z*z)
    theta = np.arctan2(np.sqrt(x*x + y*y),z)
    phi = np.arctan2(y,x)
    return r,theta,phi
B_means = []
B_stds = []
Bxy_means = []
Bxy_stds = []
zenith_means = []
zenith_stds = []
azimuth_means = []
azimuth_stds = []
for irot in rotations[:]:
    magnetometer_reading = df[[irot]].values    
    # print("magnetometer reading",magnetometer_reading)
    angles = []
    zenith = []
    azimuth = []
    B = []
    Bxy = []
    for imeas in magnetometer_reading[:]:
        nx,ny,nz = [float(jx.strip("[").strip("]")) for ix in imeas for jx in ix.split(",")]
        # print(tilt_angle(nx,ny,nz)*180/np.pi)
        angles.append(tilt_angle(nx,ny,nz)*180/np.pi)
        zenith.append(to_spherical(nx,ny,nz)[1]*180/np.pi)
        azimuth.append(to_spherical(nx,ny,nz)[2]*180/np.pi)
        B.append(to_spherical(nx,ny,nz)[0]*10**6)
        Bxy.append(np.sqrt(nx**2+ny**2))
    # print(f"Rotation {irot} mean {np.mean(angles):.2f} std {np.std(angles):.2f}")s
    # print(f"Rotation {irot} mean zen {np.mean(zenith):.2f} std {np.std(zenith):.2f}"+
    #       f" mean azi {np.mean(azimuth):.2f} std {np.std(azimuth):.2f}"+
    #       f" mean B {np.mean(B):.2f} std {np.std(B):.2f}")
    azimuth = [i+360 if i<0 else i for i in azimuth]
    azimuth_means.append(np.mean(azimuth))
    azimuth_stds.append(np.std(azimuth))
    zenith_means.append(np.mean(zenith))
    zenith_stds.append(np.std(azimuth))
    B_means.append(np.mean(B))
    B_stds.append(np.std(B))
    Bxy_means.append(np.mean(Bxy))
    Bxy_stds.append(np.std(Bxy))

# print(df.iloc[:,:1].values)
def plot_magnetometer_variation(B,B_stds,zenith,zenith_stds,azimuth,azimuth_stds):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    ax.errorbar(rotations,azimuth,yerr=azimuth_stds,fmt="o",label="azimuth")
    # ax.errorbar(rotations,B,yerr=B_stds,fmt="o",label="B")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
    ax.set_ylabel(r" azimuth of B [$^{\circ}$]", fontsize=22)
    ax.grid(True,alpha=0.6)
    ax.set_xticks(np.linspace(0,360,9))
    ax.set_yticks(np.linspace(0,360,9))
    ax.set_aspect('equal')
    ax.legend()
    ax.set_xlim(0,360)
    ax.set_ylim(0,360)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/magnetometerxyz_fine.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometerxyz_fine.pdf",transparent=False,bbox_inches='tight')
    plt.close()
plot_magnetometer_variation(B_means,B_stds,zenith_means,zenith_stds,azimuth_means,azimuth_stds)



def plot_Bxy(Bxy_means, Bx_stds):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    ax.errorbar(rotations,Bxy_means,yerr=Bxy_stds,fmt="o",label="azimuth")
    # ax.errorbar(rotations,B,yerr=B_stds,fmt="o",label="B")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
    ax.set_ylabel(r" Bxy[$\mu$T$]", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    # ax.set_aspect('equal')
    ax.legend()
    # ax.set_xlim(0,360)
    # ax.set_ylim(0,360)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/magnetometerBxy_fine.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometerBxy_fine.pdf",transparent=False,bbox_inches='tight')
    plt.close()
plot_Bxy(Bxy_means,Bxy_stds)