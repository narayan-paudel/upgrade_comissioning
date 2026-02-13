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
df_mDOM_Utah = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/utah/mDOM_North_UtahScan_PMT8_16North.txt",header=0,sep=" ")
df_mDOM_117_2 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250908XNorth_Lab117.txt",header=0,sep=" ")



# print(df_mDOMMB_yNorth.head())

def get_mean_B(df):
    '''
    returns field in muT and degrees
    '''
    bx,by,bz = df[["bx","by","bz"]].values.T
    r,theta,phi = to_spherical_list(bx,by,bz)
    if np.rad2deg(np.mean(phi))<0:
        print(f"negative  phi detected {np.rad2deg(np.mean(phi))}")
    return np.mean(bx)*10**6,np.std(bx)*10**6,np.mean(by)*10**6,np.std(by)*10**6,np.mean(bz)*10**6,np.std(bz)*10**6,\
        np.mean(r)*10**6,np.std(r)*10**6,np.rad2deg(np.mean(theta)),np.rad2deg(np.std(theta)),np.rad2deg(np.mean(phi)),np.rad2deg(np.std(phi))


def get_mean_B_rolling(df,angle):
    '''
    returns field in muT and degrees
    '''
    bx,by,bz = df[["bx","by","bz"]].values[angle:angle+10].T
    # print(f"bx {bx} by {by} bz {bz}")
    angle_list = []
    r,theta,phi = to_spherical_list(bx,by,bz)
    return np.mean(bx)*10**6,np.std(bx)*10**6,np.mean(by)*10**6,np.std(by)*10**6,np.mean(bz)*10**6,np.std(bz)*10**6,\
        np.mean(r)*10**6,np.std(r)*10**6,np.rad2deg(np.mean(theta)),np.rad2deg(np.std(theta)),np.rad2deg(np.mean(phi)),np.rad2deg(np.std(phi))

# get_mean_B_rolling(df_mDOM_rooftop,10)

angles = [i*10 for i in range(0,37)]

# for iangle in angles:
#     print(iangle)
#     print(df_mDOM_rooftop[["bx","by","bz"]].values[iangle:iangle+10])
# print(angles)


def get_mean_g(df):
    '''
    returns field in ms^-2 and degrees
    '''
    gx,gy,gz = df[["gx","gy","gz"]].values.T
    r,theta,phi = to_spherical_list(gx,gy,gz)
    return np.mean(gx),np.std(gx),np.mean(gy),np.std(gy),np.mean(gz),np.std(gz),np.mean(r),np.std(r),np.rad2deg(np.mean(theta)),\
        np.rad2deg(np.std(theta)),np.rad2deg(np.mean(phi)),np.rad2deg(np.std(phi))
# def fit_circle(x, y):
#     """
#     Fit a circle (x - x0)^2 + (y - y0)^2 = r^2
#     using least squares.
#     Returns (x0, y0, r).
#     """
#     A = np.c_[2*x, 2*y, np.ones_like(x)]
#     b = x**2 + y**2
#     sol, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
#     x0, y0, c = sol
#     r = np.sqrt(c + x0**2 + y0**2)
#     return x0, y0, r

import numpy as np
from scipy import optimize

def calc_R(x, y, xc, yc):
    """
    Calculate the distance of each point from the center (xc, yc).
    """
    return np.sqrt((x - xc)**2 + (y - yc)**2)

def residuals(params, x, y):
    xc, yc, r = params
    return calc_R(x, y, xc, yc) - r

def fit_circle(x, y):
    """
    Circle fit using linear least squares.
    (x - x0)^2 + (y - y0)^2 = r^2
    Returns: (x0, y0, r)
    """
    x = np.asarray(x)
    y = np.asarray(y)
    #center estimate
    x_m = np.mean(x)
    y_m = np.mean(y)
    #least square approximation
    center, ier = optimize.leastsq(residuals, (x_m, y_m, np.std(x)), args=(x, y))
    xc, yc, r = center
    return xc, yc, r

def calibrate_plane(x_list,y_list):
    """
    x_list: list of X magnetometer samples
    y_list: list of Y magnetometer samples
    Returns: calibrated [X,Y] data.
    """
    x, y = x_list, y_list

    # Fit circle to XY data
    x0, y0, r = fit_circle(x, y)
    print(x0, y0, r)
    
    # Offsets (hard-iron correction)
    offsets = np.array([x0, y0])
    
    # Scale factors (soft-iron correction → average radius)
    scale = r / np.std(np.sqrt((x - x0)**2 + (y - y0)**2))
    
    # Apply correction
    x_corr = (x - x0) * scale
    y_corr = (y - y0) * scale

    return x_corr, y_corr, offsets, scale

def plot_xy_360(df_list,df_labels,MB):
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
    ax.set_ylim(-80,50)
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

# plot_orientation_360([df_mDOM_rooftop,df_mDOM_rooftop_pos2,df_mDOM_rooftop_away],["rooftop_near_telescope","rooftop_near_telescope_pos2","rooftop_away_telescope",],MB="mDOM")
# plot_orientation_360([df_mDOM_317,df_mMB_317],["mDOM MB 317", "mMB 317"],MB="mDOM_room")
# plot_xy_360([df_mDOM_117],["mDOM MB 117"],MB="mDOM_mb_117")
# plot_xy_360([df_mDOM_117,df_mDOM_117_repeat,df_mDOM_117_pos2,df_mDOM_117_pos3,df_mDOM_117_floor],["mDOM MB 117","mDOM MB 117 Repeat","mDOM MB 117 Pos2","mDOM MB 117 Pos3","mDOM MB 117 Floor"],MB="mDOM_lab")
# plot_xy_360([df_mDOM_117],["mDOM MB 117"],MB="mDOM_mb_117")
# plot_xy_360([df_mDOM_rooftop,df_mDOM_rooftop_pos2,df_mDOM_rooftop_away,df_mMB_rooftop],["rooftop_near_telescope","rooftop_near_telescope_pos2","rooftop_away_telescope","rooftop_near_telescope_mMB"],MB="mDOM_rooftop")
# plot_xy_360([df_mDOM_rooftop_2,df_mDOM_rooftop_away_2],["rooftop_near_telescope","rooftop_away_telescope"],MB="mDOM_rooftop_sep8")
# plot_xy_360([df_mDOM_rooftop_away],["rooftop_away_telescope",],MB="mDOM_rooftop")
plot_xy_360([df_mDOM_Utah],["utah",],MB="mDOM_utah")
plot_xy_360([df_mDOM_117_2],["117"],MB="mb_117")




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
B_Utah = 50.77 #muT
inclination_Utah = 90 + 65.57 #degrees Down
declination_Utah = 90 - (10.75) #degrees -West (+ve would be east)

B_gallallee = 48.1428 #muT
inclination_gallallee = 90 + 61.5043 #degrees Down
declination_gallallee = 90 - (-3.4068) #degrees -West (+ve would be east)


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
        ax.axhline(reference_value,0,1,ls="--",lw=2.5,label=f"B$_{{geo}}$ ({reference_value:.1f} ${{\mu}}$T)",alpha=1.0)
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

plot_r_360([df_mDOM_Utah],["utah",],MB="mDOM_utah",reference_value=B_Utah)
plot_r_360([df_mDOM_117_2],["117"],MB="mb_117",reference_value=B_gallallee)

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
        ax.axhline(reference_value,0,1,ls="--",lw=2.5,label=f"$\delta$ ({reference_value:.1f}\u00b0)",alpha=1.0)
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

plot_theta_360([df_mDOM_Utah],["utah",],MB="mDOM_utah",reference_value=inclination_Utah)
plot_theta_360([df_mDOM_117_2],["117"],MB="mb_117",reference_value=inclination_gallallee)

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

plot_phi_360([df_mDOM_Utah],["utah",],MB="mDOM_utah")
plot_phi_360([df_mDOM_117_2],["117"],MB="mb_117")












bx = df_mDOM_Utah["bx"].values
by = df_mDOM_Utah["by"].values
bz = df_mDOM_Utah["bz"].values

b_map = {}

b_map["bx_step1"] = bx[:50] #PMT 8 North
b_map["by_step1"] = by[:50] #PMT 8 North
b_map["bz_step1"] = bz[:50] #PMT 8 North 
b_map["bx_step2"] = bx[60:110] #PMT 9 North
b_map["by_step2"] = by[60:110] #PMT 9 North
b_map["bz_step2"] = bz[60:110] #PMT 9 North
b_map["bx_step3"] = bx[150:200] #PMT 10 North
b_map["by_step3"] = by[150:200] #PMT 10 North
b_map["bz_step3"] = bz[150:200] #PMT 10 North
b_map["bx_step4"] = bx[290:340] #PMT 11 North
b_map["by_step4"] = by[290:340] #PMT 11 North
b_map["bz_step4"] = bz[290:340] #PMT 11 North
b_map["bx_step5"] = bx[380:430] #PMT 4 North
b_map["by_step5"] = by[380:430] #PMT 4 North
b_map["bz_step5"] = bz[380:430] #PMT 4 North
b_map["bx_step6"] = bx[500:550] #PMT 5 North
b_map["by_step6"] = by[500:550] #PMT 5 North
b_map["bz_step6"] = bz[500:550] #PMT 5 North
b_map["bx_step7"] = bx[650:700] #PMT 6 North
b_map["by_step7"] = by[650:700] #PMT 6 North
b_map["bz_step7"] = bz[650:700] #PMT 6 North
b_map["bx_step8"] = bx[750:800] #PMT 7 North
b_map["by_step8"] = by[750:800] #PMT 7 North
b_map["bz_step8"] = bz[750:800] #PMT 7 North
b_map["bx_step9"] = bx[860:910] #PMT 8 North
b_map["by_step9"] = by[860:910] #PMT 8 North
b_map["bz_step9"] = bz[860:910] #PMT 8 North
b_map["bx_step10"] = bx[1000:1050] #PMT 9 North
b_map["by_step10"] = by[1000:1050] #PMT 9 North
b_map["bz_step10"] = bz[1000:1050] #PMT 9 North
b_map["bx_step11"] = bx[1100:1150] #PMT 10 North
b_map["by_step11"] = by[1100:1150] #PMT 10 North
b_map["bz_step11"] = bz[1100:1150] #PMT 10 North
b_map["bx_step12"] = bx[1200:1250] #PMT 11 North
b_map["by_step12"] = by[1200:1250] #PMT 11 North
b_map["bz_step12"] = bz[1200:1250] #PMT 11 North

bx_mean_list = []
by_mean_list = []
bz_mean_list = []

# step_list = [1,2,3,4,5,6,7,8,9,10,11,12,13]
# step_list = [1,2,3,4,5,6,7,8,9,10,11,12]
# step_list = [2,3,4,5,6,7,8,9,10]
# step_list = [2,3,4,5,6,7,8,9,10] #first and last repeat to close circle
step_list = [2,3,4,5,6,7,8,9] #first and last doesnot repeat to close circle Start from PMT 9
# step_list = [1,2,3,4,5,6,7,8] #first and last doesnot repeat to close circle start from PMT 8

for istep in step_list:
    bx_mean = np.mean(b_map[f"bx_step{istep}"])
    by_mean = np.mean(b_map[f"by_step{istep}"])
    bz_mean = np.mean(b_map[f"bz_step{istep}"])
    bx_mean_list.append(bx_mean)
    by_mean_list.append(by_mean)
    bz_mean_list.append(bz_mean)

print(bx_mean_list)
print(by_mean_list)
print(bz_mean_list)


mdom = "mDOM_Utah"
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


plot_xyz(bx_mean_list[:], by_mean_list[:], bz_mean_list[:], label="coarse")
# plot_xy_calibrated(bx[:], by[:],label="fine")
# plot_xy_calibrated(bx[100:1050], by[100:1050],label="fine")#first and last doesnot repeat to close circle Start from PMT 8
# plot_xy_calibrated(bx[100:900], by[100:900],label="fine")#first and last doesnot repeat to close circle Start from PMT 9
# plot_xy_calibrated(bx[:900], by[:900],label="fine")#first and last doesnot repeat to close circle Start from PMT 8


def plot_corrected_heading_360(bx_mean_list,by_mean_list,label=""):
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
    ncolor = 0
    bx_calibrated, by_calibrated = corrected_ellipse(bx_mean_list, by_mean_list)
    headings_original = [np.rad2deg(np.arctan2(by, bx)) for bx, by in zip(bx_mean_list, by_mean_list)]
    print(f"headings original: {headings_original}")
    headings_original = [i+360 if i<-1.0 else i for i in headings_original]
    print(f"headings original: {headings_original}")
    headings_corrected = [np.rad2deg(np.arctan2(byc, bxc)) for bxc,byc in zip(bx_calibrated, by_calibrated)]
    print(f"headings corrected: {headings_corrected}")
    headings_corrected = [i+360 if i<-1.0 else i for i in headings_corrected]
    # headings_corrected = [i+360 if 0<i<42 else i for i in headings_corrected]
    # headings_corrected = [i+360 if 0<i<34.0 else i for i in headings_corrected]
    # print(f"headings original: {headings_original}")
    print(f"headings corrected: {headings_corrected}")
    ###################################
    if label == "coarse":
        ax.plot(((np.asarray(step_list))-1)*45, headings_original, "-o", c=colorsCustom[ncolor], label=f"{'raw'}", alpha=1)
        ax.plot(((np.asarray(step_list))-1)*45, headings_corrected, "--o", c=colorsCustom[ncolor+2], label=f"{'calibrated'}", alpha=1)
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
    # ax.set_xlim(0,360)
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
    # ax.set_aspect('equal')
    ax.legend(ncols=2,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" mDOM rotation [$^{\circ}$]", fontsize=20)
    ax.set_ylabel(r" $\phi$ [$^{\circ}$]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_Bxy_heading_Utah_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_Bxy_heading_Utah_{label}.pdf",transparent=False,bbox_inches='tight')
    plt.close()
plot_corrected_heading_360(bx_mean_list, by_mean_list,label="coarse")
# plot_corrected_heading_360(bx, by,label="fine")
