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

#measurement taken on August 20, 2025 for four orientation of magnetometer y: North, West, South, East
sensor_readings_mDOM = sorted(glob.glob("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250820Y*.txt"))
sensor_readings_mDOM_yNorth = sorted(glob.glob("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250820YN*.txt"))
sensor_readings_mDOM_yWest = sorted(glob.glob("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250820YW*.txt"))
sensor_readings_mDOM_ySouth = sorted(glob.glob("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250820YS*.txt"))
sensor_readings_mDOM_yEast = sorted(glob.glob("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250820YE*.txt"))
sensor_readings_mMB = glob.glob("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmMB_20250820Y*.txt")
sensor_readings_mMB_yNorth = sorted(glob.glob("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmMB_20250820YN*.txt"))
sensor_readings_mMB_yWest = sorted(glob.glob("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmMB_20250820YW*.txt"))
sensor_readings_mMB_ySouth = sorted(glob.glob("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmMB_20250820YS*.txt"))
sensor_readings_mMB_yEast = sorted(glob.glob("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmMB_20250820YE*.txt"))

which_run = 0 #0 or 1
print(sensor_readings_mDOM_yNorth[which_run])
print(sensor_readings_mDOM_yWest[which_run])
print(sensor_readings_mDOM_ySouth[which_run])
print(sensor_readings_mDOM_yEast[which_run])



import pandas as pd

df_mDOMMB_yNorth = pd.read_csv(sensor_readings_mDOM_yNorth[which_run],header=0,sep=" ")
df_mDOMMB_yWest = pd.read_csv(sensor_readings_mDOM_yWest[which_run],header=0,sep=" ")
df_mDOMMB_ySouth = pd.read_csv(sensor_readings_mDOM_ySouth[which_run],header=0,sep=" ")
df_mDOMMB_yEast = pd.read_csv(sensor_readings_mDOM_yEast[which_run],header=0,sep=" ")

df_mDOMMB_yNorthUp = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250820YNorth3.txt",header=0,sep=" ")
df_mDOMMB_yNorthDown = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250820YNorth3Invert.txt",header=0,sep=" ")


df_mMB_yNorth = pd.read_csv(sensor_readings_mMB_yNorth[which_run],header=0,sep=" ")
df_mMB_yWest = pd.read_csv(sensor_readings_mMB_yWest[which_run],header=0,sep=" ")
df_mMB_ySouth = pd.read_csv(sensor_readings_mMB_ySouth[which_run],header=0,sep=" ")
df_mMB_yEast = pd.read_csv(sensor_readings_mMB_yEast[which_run],header=0,sep=" ")

df_mMB_yNorthUp = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmMB_20250820YNorth3.txt",header=0,sep=" ")
df_mMB_yNorthDown = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmMB_20250820YNorth3Invert.txt",header=0,sep=" ")





# print(df_mDOMMB_yNorth.head())

def get_mean_B(df):
    '''
    returns field in muT and degrees
    '''
    bx,by,bz = df[["bx","by","bz"]].values.T
    r,theta,phi = to_spherical_list(bx,by,bz)
    return np.mean(bx)*10**6,np.std(bx)*10**6,np.mean(by)*10**6,np.std(by)*10**6,np.mean(bz)*10**6,np.std(bz)*10**6,\
        np.mean(r)*10**6,np.std(r)*10**6,np.rad2deg(np.mean(theta)),np.rad2deg(np.std(theta)),np.rad2deg(np.mean(phi)),np.rad2deg(np.std(phi))


def get_mean_g(df):
    '''
    returns field in ms^-2 and degrees
    '''
    gx,gy,gz = df[["gx","gy","gz"]].values.T
    r,theta,phi = to_spherical_list(gx,gy,gz)
    return np.mean(gx),np.std(gx),np.mean(gy),np.std(gy),np.mean(gz),np.std(gz),np.mean(r),np.std(r),np.rad2deg(np.mean(theta)),\
        np.rad2deg(np.std(theta)),np.rad2deg(np.mean(phi)),np.rad2deg(np.std(phi))

def plot_orientation(df_list,dir_list,MB):
    fig = plt.figure(figsize=(8,5*3))
    gs = gridspec.GridSpec(nrows=3,ncols=1, figure=fig)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[2,0])
    for df,dir in zip(df_list,dir_list):
        bx,bx_std,by,by_std,bz,bz_std,r,r_std,theta,theta_std,phi,phi_std = get_mean_B(df)
        print(f"{dir} Bx {bx:.2f}+-{bx_std:.2f} By {by:.2f}+-{by_std:.2f} Bz {bz:.2f}+-{bz_std:.2f} B {r:.2f}+-{r_std:.2f} theta {theta:.2f}+-{theta_std:.2f} phi {phi:.2f}+-{phi_std:.2f}")
        ax1.errorbar(dir,r,yerr=r_std,fmt="o",label=f"{""}")
        ax2.errorbar(dir,theta,yerr=theta_std,fmt="o",label=f"{""}")
        ax3.errorbar(dir,phi,yerr=phi_std,fmt="o",label=f"{""}")
    ax1.axhline(B_gallallee,0,1,ls="--",lw=2.5,label=f"B$_{{geo}}$ ({B_gallallee:.1f} ${{\mu}}$T)",alpha=1.0)
    ax2.axhline(inclination_gallallee,0,1,ls="--",lw=2.5,label=f"$\delta$ ({inclination_gallallee:.1f}\u00b0)",alpha=1.0)
    ax3.axhline(declination_gallallee,0,1,ls="--",lw=2.5,label=f"$I$ ({declination_gallallee:.1f}\u00b0)",alpha=1.0)
    for ax in [ax1,ax2,ax3]:
        ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
        ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)        
        ax.grid(True,alpha=0.6)
        # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
        # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
        # ax.set_xlim(0,360)
        # ax.set_xticks(np.linspace(0,360,9))
        # ax.set_yticks(np.linspace(0,360,9))
        # ax.set_aspect('equal')
        ax.legend(ncols=5)
    
    ax2.set_yticks(np.linspace(0,360,37))
    ax3.set_yticks(np.linspace(0,360,19))
    ax1.set_ylabel(r" B [$\mu$T]", fontsize=22)
    ax2.set_ylabel(r" $\theta$ [$^{\circ}$]", fontsize=22)
    ax3.set_ylabel(r" $\phi$ [$^{\circ}$]", fontsize=22)
    # ax1.set_ylim(35,130)
    ax1.set_ylim(0,130)
    # ax3.set_yticks(np.linspace(0,360,5))
    ax2.set_ylim(0,180)
    # ax2.set_ylim(-1.5,1.5)
    # ax3.set_ylim(8,11)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/orientation_with_B{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_B{MB}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_orientation(df_list=[df_mDOMMB_yNorth,df_mDOMMB_yWest,df_mDOMMB_ySouth,df_mDOMMB_yEast],dir_list=["North","West","South","East"],MB="mDOMMB")
# plot_orientation(df_list=[df_mMB_yNorth,df_mMB_yWest,df_mMB_ySouth,df_mMB_yEast],dir_list=["North","West","South","East"],MB="mMB")
# plot_orientation(df_list=[df_mDOMMB_yNorthUp,df_mDOMMB_yNorthDown],dir_list=["Upright","Inverted"],MB="InvertmDOM")
# plot_orientation(df_list=[df_mMB_yNorthUp,df_mMB_yNorthDown],dir_list=["Upright","Inverted"],MB="InvertmMB")


def plot_orientation_g(df_list,dir_list,MB):
    fig = plt.figure(figsize=(8,5*3))
    gs = gridspec.GridSpec(nrows=3,ncols=1, figure=fig)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[2,0])
    for df,dir in zip(df_list,dir_list):
        gx,gx_std,gy,gy_std,gz,gz_std,r,r_std,theta,theta_std,phi,phi_std = get_mean_g(df)
        print(f"{dir} gx {gx:.2f}+-{gx_std:.2f} gy {gy:.2f}+-{gy_std:.2f} gz {gz:.2f}+-{gz_std:.2f} B {r:.2f}+-{r_std:.2f} theta {theta:.2f}+-{theta_std:.2f} phi {phi:.2f}+-{phi_std:.2f}")
        ax1.errorbar(dir,r,yerr=r_std,fmt="o",label=f"{""}")
        ax2.errorbar(dir,theta,yerr=theta_std,fmt="o",label=f"{""}")
        ax3.errorbar(dir,phi,yerr=phi_std,fmt="o",label=f"{""}")
    # ax1.axhline(B_gallallee,0,1,ls="--",lw=2.5,label=f"B$_{{geo}}$ ({B_gallallee:.1f} ${{\mu}}$T)",alpha=1.0)
    # ax2.axhline(inclination_gallallee,0,1,ls="--",lw=2.5,label=f"$\delta$ ({inclination_gallallee:.1f}\u00b0)",alpha=1.0)
    # ax3.axhline(declination_gallallee,0,1,ls="--",lw=2.5,label=f"$I$ ({declination_gallallee:.1f}\u00b0)",alpha=1.0)
    for ax in [ax1,ax2,ax3]:
        ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
        ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)        
        ax.grid(True,alpha=0.6)
        # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
        # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
        # ax.set_xlim(0,360)
        # ax.set_xticks(np.linspace(0,360,9))
        # ax.set_yticks(np.linspace(0,360,9))
        # ax.set_aspect('equal')
        ax.legend(ncols=5)
    
    ax2.set_yticks(np.linspace(0,360,37))
    ax3.set_yticks(np.linspace(0,360,19))
    ax1.set_ylabel(r" g [ms$^{-2}$]", fontsize=22)
    ax2.set_ylabel(r" $\theta$ [$^{\circ}$]", fontsize=22)
    ax3.set_ylabel(r" $\phi$ [$^{\circ}$]", fontsize=22)
    # ax1.set_ylim(35,130)
    ax1.set_ylim(0,12)
    # ax3.set_yticks(np.linspace(0,360,5))
    ax2.set_ylim(0,180)
    # ax2.set_ylim(-1.5,1.5)
    # ax3.set_ylim(8,11)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/orientation_with_g{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_g{MB}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_orientation_g(df_list=[df_mDOMMB_yNorth,df_mDOMMB_yWest,df_mDOMMB_ySouth,df_mDOMMB_yEast],dir_list=["North","West","South","East"],MB="mDOMMB")
plot_orientation_g(df_list=[df_mMB_yNorth,df_mMB_yWest,df_mMB_ySouth,df_mMB_yEast],dir_list=["North","West","South","East"],MB="mMB")
plot_orientation_g(df_list=[df_mDOMMB_yNorthUp,df_mDOMMB_yNorthDown],dir_list=["Upright","Inverted"],MB="InvertmDOM")
plot_orientation_g(df_list=[df_mMB_yNorthUp,df_mMB_yNorthDown],dir_list=["Upright","Inverted"],MB="InvertmMB")




