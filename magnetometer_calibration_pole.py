#!/usr/bin/env python

import os
import glob

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pathlib import Path
home = str(Path.home())

from utils import to_spherical, tilt_angle, get_means, get_sub_df, to_spherical_list, to_spherical_list_rotated,meas_from_df
from calibrate_magnetometer_2d import corrected_ellipse

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = home+"/research_ua/icecube/upgrade/upgrade_comissioning/plots/"
# print(degg_list)

colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
B_gallallee = 48.1428 #muT
inclination_gallallee = 90 + 61.5043 #degrees Down
declination_gallallee = 90 - (-3.4068) #degrees -West (+ve would be east)


drive = "/Users/epaudel/Library/CloudStorage/OneDrive-TheUniversityofAlabama/"
accelerometer_readings = drive+"research_notes/references/upgrade_comissioning/mDOM_orientations/accelerometer_measurement.xlsx"

#measurement taken on August 25, 2025 for magnetometer yNorth --> 360 rotation
import pandas as pd

#measurement at mDOM at Utah on Oct 09 2025, PMT 8/16 facing North (the water tank), rotating anticlockwise when looking from top
# df_mDOM_Utah = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/utah/mDOM_North_UtahScan.txt",header=0,sep=" ")
df_mDOM_pole_run1 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/pole_rotation/sensorsmDOMMB_20251215_run1.txt",header=0,sep=" ")
df_mDOM_pole_run2 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/pole_rotation/sensorsmDOMMB_20251215_run2.txt",header=0,sep=" ")
df_mDOM_pole_run3 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/pole_rotation/sensorsmDOMMB_20251215_run3.txt",header=0,sep=" ")
df_mDOM_pole_run4 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/pole_rotation/sensorsmDOMMB_20251215_run4.txt",header=0,sep=" ")


run = 3
df_mDOM = df_mDOM_pole_run1

def to_spherical(x,y,z):
    '''
    converts cartesian coordinates (x,y,z) to spherical (r,theta,phi)
    '''
    r = np.sqrt(x*x+y*y+z*z)
    theta = np.arctan2(np.sqrt(x*x + y*y),z)
    phi = np.arctan2(y,x)
    return r,theta,phi

def spherical_lists(x_list,y_list,z_list):
    r_list = []
    θ_list = []
    φ_list = []
    for x,y,z in zip(x_list,y_list,z_list):
        r,theta,phi = to_spherical(x,y,z)
        r_list.append(r)
        θ_list.append(theta)
        φ_list.append(phi)
    return r_list,θ_list,φ_list
B_pole = 50.77 #muT
inclination_pole = 90 + 65.57 #degrees Down
declination_pole = 90 - (10.75) #degrees -West (+ve would be east)

# B_gallallee = 48.1428 #muT
# inclination_gallallee = 90 + 61.5043 #degrees Down
# declination_gallallee = 90 - (-3.4068) #degrees -West (+ve would be east)

def plot_xy_360(df_list,df_labels,MB):
    # fig = plt.figure(figsize=(8,5))
    fig = plt.figure(figsize=(8*10,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
    ncolor = 0
    for df,df_label in zip(df_list,df_labels):
        bx = df["bx"].values *10**6  #convert to microTesla
        # bx = bx[:60] #convert to microTesla
        by = df["by"].values *10**6  #convert to microTesla
        # by = by[:60] #convert to microTesla
        bz = df["bz"].values *10**6  #convert to microTesla
        # bz = bz[:60] #convert to microTesla
        ax.plot(range(len(bx)), bx,label=f"{'B$_{x}$'}", alpha=1)
        ax.plot(range(len(by)), by,label=f"{'B$_{y}$'}", alpha=1)
        ax.plot(range(len(bz)), bz,label=f"{'B$_{z}$'}", alpha=1)
        # ax.plot(bx_calibrated, by_calibrated, "--o",c=colorsCustom[ncolor+2], label=f"{df_label} (calibrated)", alpha=1)
        ncolor+=1
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xlim(0,360)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    # ax.set_aspect('equal')
    ax.legend(loc=(0.35,0.25),ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" # of samples", fontsize=20)
    ax.set_ylabel(r" $B_i$ [$\mu$T]", fontsize=20)
    # ax.set_ylim(-80,50)
    # ax.set_xlim(1100,1400)
    # ax.set_ylim(-55,55)
    # ax.set_xlim(-55,55)
    # ax3.set_yticks(np.linspace(0,360,5))
    # ax.set_ylim(0,180)
    # ax2.set_ylim(-1.5,1.5)
    # ax3.set_ylim(8,11)
    # ax.set_xlim(1300,1500)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/orientation_with_B_{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_B_{MB}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_xy_360([df_mDOM],["pole",],MB="mDOM_pole")




def plot_r_360(df_list,df_labels,MB,reference_value):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
    ncolor = 0
    for df,df_label in zip(df_list,df_labels):
        bx = df["bx"].values *10**6  #convert to microTesla
        # bx = bx[:60] #convert to microTesla
        by = df["by"].values *10**6  #convert to microTesla
        # by = by[:60] #convert to microTesla
        bz = df["bz"].values *10**6  #convert to microTesla
        r,theta,phi = spherical_lists(bx,by,bz)

        # bz = bz[:60] #convert to microTesla
        ax.plot(range(len(r)), r,label=f"{'B'}", alpha=1)
        # ax.axhline(reference_value,0,1,ls="--",lw=2.5,label=f"B$_{{geo}}$ ({reference_value:.1f} ${{'\u03bc'}}$T)",alpha=1.0)
        # ax.plot(bx_calibrated, by_calibrated, "--o",c=colorsCustom[ncolor+2], label=f"{df_label} (calibrated)", alpha=1)
        ncolor+=1    
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.legend(loc=(0.35,0.25),ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" # of samples", fontsize=20)
    ax.set_ylabel(r" $B_i$ [$\mu$T]", fontsize=20)
    ax.set_ylim(0,90)
    plt.savefig(plotFolder+f"/orientation_with_B_r_{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_B_r_{MB}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_r_360([df_mDOM],["pole",],MB="mDOM_pole",reference_value=B_pole)

def plot_theta_360(df_list,df_labels,MB,reference_value):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
    ncolor = 0
    for df,df_label in zip(df_list,df_labels):
        bx = df["bx"].values *10**6  #convert to microTesla
        # bx = bx[:60] #convert to microTesla
        by = df["by"].values *10**6  #convert to microTesla
        # by = by[:60] #convert to microTesla
        bz = df["bz"].values *10**6  #convert to microTesla
        r,theta,phi = spherical_lists(bx,by,bz)
        # bz = bz[:60] #convert to microTesla
        ax.plot(range(len(theta)), np.rad2deg(np.asarray(theta)), "-o", c=colorsCustom[0], label=f"$\\theta$\u00b0", alpha=1)
        # ax.axhline(reference_value,0,1,ls="--",lw=2.5,label=f"$\\delta$ ({reference_value:.1f}\u00b0)",alpha=1.0)
        # ax.plot(bx_calibrated, by_calibrated, "--o",c=colorsCustom[ncolor+2], label=f"{df_label} (calibrated)", alpha=1)
        ncolor+=1    
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.legend(loc=(0.35,0.25),ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_ylim(0,180)
    ax.set_xlabel(r" # of samples", fontsize=20)
    ax.set_ylabel(f" $\\theta$\u00b0", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_B_Theta_{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_B_Theta_{MB}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_theta_360([df_mDOM],["pole",],MB="mDOM_pole",reference_value=inclination_pole)

def plot_phi_360(df_list,df_labels,MB):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
    ncolor = 0
    for df,df_label in zip(df_list,df_labels):
        bx = df["bx"].values *10**6  #convert to microTesla
        # bx = bx[:60] #convert to microTesla
        by = df["by"].values *10**6  #convert to microTesla
        # by = by[:60] #convert to microTesla
        bz = df["bz"].values *10**6  #convert to microTesla
        r,theta,phi = spherical_lists(bx,by,bz)
        # bz = bz[:60] #convert to microTesla
        ax.plot(range(len(phi)), np.rad2deg(np.asarray(phi)), "-o", c=colorsCustom[0], label=f"$\\phi$\u00b0", alpha=1)
        # ax.axhline(inclination_Utah,0,1,ls="--",lw=2.5,label=f"$\delta$ ({inclination_Utah:.1f}\u00b0)",alpha=1.0)
        # ax.plot(bx_calibrated, by_calibrated, "--o",c=colorsCustom[ncolor+2], label=f"{df_label} (calibrated)", alpha=1)
        ncolor+=1    
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.legend(loc=(0.35,0.25),ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_ylim(-180,180)
    ax.set_xlabel(r" # of samples", fontsize=20)
    ax.set_ylabel(f" $\\phi$\u00b0", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_B_Phi_{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_B_Phi_{MB}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_phi_360([df_mDOM],["pole",],MB="mDOM_pole")


bx = df_mDOM["bx"].values
by = df_mDOM["by"].values
bz = df_mDOM["bz"].values

b_map = {}
mdom = ""

def plot_xyz(mag_x, mag_y, mag_z, label=""):
    mag_x = np.array(mag_x)*10**6 #convert to microTesla
    mag_y = np.array(mag_y)*10**6 #convert to microTesla
    mag_z = np.array(mag_z)*10**6 #convert to microTesla
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0],projection='3d')
    ax.scatter3D(mag_x, mag_y,mag_z, "-o", c="b", label=f"{"raw"}", alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.set_aspect('equal')
    ax.legend(loc="lower left",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" $B_x$ [$\mu$T]", fontsize=16)
    ax.set_ylabel(r" $B_y$ [$\mu$T]", fontsize=16)
    ax.set_zlabel(r" $B_z$ [$\mu$T]", fontsize=16)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xyz_calibrated_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xyz_calibrated_{label}.pdf",transparent=False,bbox_inches='tight')
    plt.show()
    # plt.close()


# plot_xyz(bx[:], by[:],bz[:], label="coarse")
# plot_xy_calibrated(bx[:], by[:],label="fine")

def plot_xy(mag_x, mag_y, label=""):
    mag_x = np.array(mag_x)*10**6 #convert to microTesla
    mag_y = np.array(mag_y)*10**6 #convert to microTesla
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    ax.plot(mag_x, mag_y, "-o", c="b", label=f"{"raw"}", alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.set_aspect('equal')
    ax.legend(loc="lower left",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlim(-50,50)
    ax.set_ylim(-50,50)
    ax.set_xlabel(r" $B_x$ [$\mu$T]", fontsize=16)
    ax.set_ylabel(r" $B_y$ [$\mu$T]", fontsize=16)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xy_calibrated_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xy_calibrated_{label}.pdf",transparent=False,bbox_inches='tight')
    # plt.show()
    plt.close()


plot_xy(bx[:], by[:], label="fine")

b_map = {}
if run == 4:
    b_map["bx_PMT8"] = bx[200:225] #PMT 8 North
    b_map["by_PMT8"] = by[200:225] #PMT 8 North
    b_map["bz_PMT8"] = bz[200:225] #PMT 8 North
    b_map["bx_PMT8p"] = bx[300:325] #PMT 8 and 9 North
    b_map["by_PMT8p"] = by[300:325] #PMT 8 and 9 North
    b_map["bz_PMT8p"] = bz[300:325] 
    b_map["bx_PMT9"] = bx[400:425] 
    b_map["by_PMT9"] = by[400:425] 
    b_map["bz_PMT9"] = bz[400:425] 
    b_map["bx_PMT9p"] = bx[500:525] 
    b_map["by_PMT9p"] = by[500:525] 
    b_map["bz_PMT9p"] = bz[500:525] 
    b_map["bx_PMT10"] = bx[575:600] 
    b_map["by_PMT10"] = by[575:600] 
    b_map["bz_PMT10"] = bz[575:600] 
    b_map["bx_PMT10p"] = bx[625:650] 
    b_map["by_PMT10p"] = by[625:650] 
    b_map["bz_PMT10p"] = bz[625:650] 
    b_map["bx_PMT11"] = bx[700:725] 
    b_map["by_PMT11"] = by[700:725] 
    b_map["bz_PMT11"] = bz[700:725] 
    b_map["bx_PMT11p"] = bx[750:775] 
    b_map["by_PMT11p"] = by[750:775] 
    b_map["bz_PMT11p"] = bz[750:775] 
    b_map["bx_PMT4"] = bx[800:825] 
    b_map["by_PMT4"] = by[800:825] 
    b_map["bz_PMT4"] = bz[800:825] 
    b_map["bx_PMT4p"] = bx[875:900] 
    b_map["by_PMT4p"] = by[875:900] 
    b_map["bz_PMT4p"] = bz[875:900] 
    b_map["bx_PMT5"] = bx[925:950] 
    b_map["by_PMT5"] = by[925:950] 
    b_map["bz_PMT5"] = bz[925:950] 
    b_map["bx_PMT5p"] = bx[975:1000] 
    b_map["by_PMT5p"] = by[975:1000] 
    b_map["bz_PMT5p"] = bz[975:1000] 
    b_map["bx_PMT6"] = bx[1050:1075] 
    b_map["by_PMT6"] = by[1050:1075] 
    b_map["bz_PMT6"] = bz[1050:1075] 
    b_map["bx_PMT6p"] = bx[1100:1125] 
    b_map["by_PMT6p"] = by[1100:1125] 
    b_map["bz_PMT6p"] = bz[1100:1125] 
    b_map["bx_PMT7"] = bx[1150:1175] 
    b_map["by_PMT7"] = by[1150:1175] 
    b_map["bz_PMT7"] = bz[1150:1175] 
    b_map["bx_PMT7p"] = bx[1225:1250] 
    b_map["by_PMT7p"] = by[1225:1250] 
    b_map["bz_PMT7p"] = bz[1225:1250] 
    b_map["bx_PMT8"] = bx[1280:1305] 
    b_map["by_PMT8"] = by[1280:1305] 
    b_map["bz_PMT8"] = bz[1280:1305] 
    print(f"bmap {b_map}")
elif run == 3:
    b_map["bx_PMT8"] = bx[0:25] #PMT 8 North
    b_map["by_PMT8"] = by[200:225] #PMT 8 North
    b_map["bz_PMT8"] = bz[200:225] #PMT 8 North
    b_map["bx_PMT8p"] = bx[300:325] #PMT 8 and 9 North
    b_map["by_PMT8p"] = by[300:325] #PMT 8 and 9 North
    b_map["bz_PMT8p"] = bz[300:325] 
    b_map["bx_PMT9"] = bx[400:425] 
    b_map["by_PMT9"] = by[400:425] 
    b_map["bz_PMT9"] = bz[400:425] 
    b_map["bx_PMT9p"] = bx[500:525] 
    b_map["by_PMT9p"] = by[500:525] 
    b_map["bz_PMT9p"] = bz[500:525] 
    b_map["bx_PMT10"] = bx[575:600] 
    b_map["by_PMT10"] = by[575:600] 
    b_map["bz_PMT10"] = bz[575:600] 
    b_map["bx_PMT10p"] = bx[625:650] 
    b_map["by_PMT10p"] = by[625:650] 
    b_map["bz_PMT10p"] = bz[625:650] 
    b_map["bx_PMT11"] = bx[700:725] 
    b_map["by_PMT11"] = by[700:725] 
    b_map["bz_PMT11"] = bz[700:725] 
    b_map["bx_PMT11p"] = bx[750:775] 
    b_map["by_PMT11p"] = by[750:775] 
    b_map["bz_PMT11p"] = bz[750:775] 
    b_map["bx_PMT4"] = bx[800:825] 
    b_map["by_PMT4"] = by[800:825] 
    b_map["bz_PMT4"] = bz[800:825] 
    b_map["bx_PMT4p"] = bx[875:900] 
    b_map["by_PMT4p"] = by[875:900] 
    b_map["bz_PMT4p"] = bz[875:900] 
    b_map["bx_PMT5"] = bx[925:950] 
    b_map["by_PMT5"] = by[925:950] 
    b_map["bz_PMT5"] = bz[925:950] 
    b_map["bx_PMT5p"] = bx[975:1000] 
    b_map["by_PMT5p"] = by[975:1000] 
    b_map["bz_PMT5p"] = bz[975:1000] 
    b_map["bx_PMT6"] = bx[1050:1075] 
    b_map["by_PMT6"] = by[1050:1075] 
    b_map["bz_PMT6"] = bz[1050:1075] 
    b_map["bx_PMT6p"] = bx[1100:1125] 
    b_map["by_PMT6p"] = by[1100:1125] 
    b_map["bz_PMT6p"] = bz[1100:1125] 
    b_map["bx_PMT7"] = bx[1150:1175] 
    b_map["by_PMT7"] = by[1150:1175] 
    b_map["bz_PMT7"] = bz[1150:1175] 
    b_map["bx_PMT7p"] = bx[1225:1250] 
    b_map["by_PMT7p"] = by[1225:1250] 
    b_map["bz_PMT7p"] = bz[1225:1250] 
    b_map["bx_PMT8"] = bx[1280:1305] 
    b_map["by_PMT8"] = by[1280:1305] 
    b_map["bz_PMT8"] = bz[1280:1305] 
    print(f"bmap {b_map}")


bx_collect = np.concatenate((b_map["bx_PMT8"],b_map["bx_PMT8p"],b_map["bx_PMT9"],b_map["bx_PMT9p"],b_map["bx_PMT10"],b_map["bx_PMT10p"],b_map["bx_PMT11"],b_map["bx_PMT11p"],b_map["bx_PMT4"],b_map["bx_PMT4p"],b_map["bx_PMT5"],b_map["bx_PMT5p"],b_map["bx_PMT6"],b_map["bx_PMT6p"],b_map["bx_PMT7"],b_map["bx_PMT7p"],b_map["bx_PMT8"]))
by_collect = np.concatenate((b_map["by_PMT8"],b_map["by_PMT8p"],b_map["by_PMT9"],b_map["by_PMT9p"],b_map["by_PMT10"],b_map["by_PMT10p"],b_map["by_PMT11"],b_map["by_PMT11p"],b_map["by_PMT4"],b_map["by_PMT4p"],b_map["by_PMT5"],b_map["by_PMT5p"],b_map["by_PMT6"],b_map["by_PMT6p"],b_map["by_PMT7"],b_map["by_PMT7p"],b_map["by_PMT8"]))
bz_collect = np.concatenate((b_map["bz_PMT8"],b_map["bz_PMT8p"],b_map["bz_PMT9"],b_map["bz_PMT9p"],b_map["bz_PMT10"],b_map["bz_PMT10p"],b_map["bz_PMT11"],b_map["bz_PMT11p"],b_map["bz_PMT4"],b_map["bz_PMT4p"],b_map["bz_PMT5"],b_map["bz_PMT5p"],b_map["bz_PMT6"],b_map["bz_PMT6p"],b_map["bz_PMT7"],b_map["bz_PMT7p"],b_map["bz_PMT8"]))

print(f"bx collect {len(bx_collect)}")



def plot_xyz_broken(mag_x, mag_y,mag_z, label=""):
    print(f"plotting broken plots")
    mag_x = np.array(mag_x)*10**6 #convert to microTesla
    mag_y = np.array(mag_y)*10**6 #convert to microTesla
    mag_z = np.array(mag_z)*10**6 #convert to microTesla
    fig = plt.figure(figsize=(24,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    print("mag x", mag_x)
    ax.plot(range(len(mag_x)), mag_x, "-o", label=f"{"x"}", alpha=1)
    ax.plot(range(len(mag_y)), mag_y, "-o", label=f"{"y"}", alpha=1)
    ax.plot(range(len(mag_z)), mag_z, "-o", label=f"{"z"}", alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.set_aspect('equal')
    ax.legend(loc="lower left",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    # ax.set_xlim(-50,50)
    # ax.set_ylim(-50,50)
    ax.set_xlabel(r"# of bins", fontsize=16)
    ax.set_ylabel(r" $B$ [$\mu$T]", fontsize=16)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xyz_broken_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xyz_broken_{label}.pdf",transparent=False,bbox_inches='tight')
    # plt.show()
    print(f"plots saved to {plotFolder+f'/orientation_with_{mdom}xyz_broken_{label}.png'}")
    plt.close()

plot_xyz_broken(bx_collect[:], by_collect[:], bz_collect[:], label="broken")


bx_mean_list = []
by_mean_list = []
bz_mean_list = []

step_list = []
for i in [8,9,10,11,4,5,6,7]:
    step_list.append(f"{i}")
    step_list.append(f"{i}p")
step_list.append(f"{8}")
pmt_magfield_angles = {}
angles_list = [i*22.5 for i in range(17)]
# print(angles_list)
# for i,j in zip(step_list,angles_list):
#     pmt_magfield_angles[i] = j

# print(f"pmt_mag_field {pmt_magfield_angles}")


for istep in step_list:
    bx_mean = np.mean(b_map[f"bx_PMT{istep}"])
    by_mean = np.mean(b_map[f"by_PMT{istep}"])
    bz_mean = np.mean(b_map[f"bz_PMT{istep}"])
    bx_mean_list.append(bx_mean)
    by_mean_list.append(by_mean)
    bz_mean_list.append(bz_mean)

# print(bx_mean_list)
# print(by_mean_list)
# print(bz_mean_list)


mdom = "mDOM_pole"
def plot_xy_calibrated(mag_x, mag_y,label=""):
    mag_x = np.array(mag_x)*10**6 #convert to microTesla
    mag_y = np.array(mag_y)*10**6 #convert to microTesla
    mag_x_cal, mag_y_cal = corrected_ellipse(mag_x, mag_y)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    ax.plot(mag_x, mag_y, "-o", c="b", label=f"{"raw"}", alpha=1)
    ax.plot(mag_x_cal, mag_y_cal, "-o", c="r", label=f"{"cal"}", alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.set_aspect('equal')
    ax.legend(loc="lower left",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" $B_x$ [$\mu$T]", fontsize=20)
    ax.set_ylabel(r" $B_y$ [$\mu$T]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xy_calibrated_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}xy_calibrated_{label}.pdf",transparent=False,bbox_inches='tight')
    plt.close()


plot_xy_calibrated(bx_mean_list[:], by_mean_list[:],label="coarse")

def plot_corrected_heading_360(bx_mean_list,by_mean_list,label=""):
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
    ncolor = 0
    bx_calibrated, by_calibrated = corrected_ellipse(bx_mean_list, by_mean_list)
    headings_original = [np.rad2deg(np.arctan2(by, bx)) for bx, by in zip(bx_mean_list, by_mean_list)]
    # print(f"headings original: {headings_original}")
    headings_original = [i+360 if i<-1.0 else i for i in headings_original]
    # print(f"headings original: {headings_original}")
    headings_corrected = [np.rad2deg(np.arctan2(byc, bxc)) for bxc,byc in zip(bx_calibrated, by_calibrated)]
    # print(f"headings corrected: {headings_corrected}")
    headings_corrected = [i+360 if i<-1.0 else i for i in headings_corrected]
    # headings_corrected = [i+360 if 0<i<42 else i for i in headings_corrected]
    # headings_corrected = [i+360 if 0<i<34.0 else i for i in headings_corrected]
    # print(f"headings original: {headings_original}")
    # print(f"headings corrected: {headings_corrected}")
    ###################################
    if label == "coarse":
        ax.plot(angles_list, headings_original, "-o", c=colorsCustom[ncolor], label=f"{'raw'}", alpha=1)
        ax.plot(angles_list, headings_corrected, "--o", c=colorsCustom[ncolor+2], label=f"{'calibrated'}", alpha=1)
        ax.set_xticks(np.linspace(0,360,9))
    else:
        ax.plot(range(len(headings_original)), headings_original, "-o", c=colorsCustom[ncolor], label=f"{'raw'}", alpha=1)
        ax.plot(range(len(headings_corrected)), headings_corrected, "--o", c=colorsCustom[ncolor+2], label=f"{'calibrated'}", alpha=1)
    #####################################
    ncolor += 1
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    ax.set_xlim(0,360)
    ax.set_ylim(0,360)
    ###############################################
    ##############################################
    # ax.set_xticks(np.linspace(0,360,9)-roll)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    # ax.set_xticks(np.linspace(0,540,13))
    # ax.set_yticks(np.linspace(-180,180,9))
    #############################################
    #############################################
    ax.set_aspect('equal')
    ax.legend(ncols=2,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" mDOM rotation [$^{\circ}$]", fontsize=20)
    ax.set_ylabel(r" $\phi$ [$^{\circ}$]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_Bxy_heading_Pole_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_Bxy_heading_Pole_{label}.pdf",transparent=False,bbox_inches='tight')
    plt.close()
plot_corrected_heading_360(bx_mean_list, by_mean_list,label="fine")
plot_corrected_heading_360(bx_mean_list, by_mean_list,label="coarse")