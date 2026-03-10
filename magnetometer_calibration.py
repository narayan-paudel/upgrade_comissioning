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

df_mDOM_rooftop = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250825YNorth_test_Rooftop.txt",header=0,sep=" ")
df_mDOM_rooftop_pos2 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250825YNorth_test_RooftopLaptopPos2.txt",header=0,sep=" ")
df_mDOM_rooftop_away = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250825YNorth_test_RooftopAwayFromTelescope.txt",header=0,sep=" ")
df_mMB_rooftop = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmMB_20250825YNorth_test_Rooftop.txt",header=0,sep=" ")

df_mDOM_317 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250825YNorth_test.txt",header=0,sep=" ")
df_mMB_317 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmMB_20250825YNorth_test.txt",header=0,sep=" ")

df_mDOM_117 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250825YNorth_test_table_Lab117.txt",header=0,sep=" ")
df_mDOM_117_repeat = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250825YNorth_test_table_Lab117Repeat.txt",header=0,sep=" ")
df_mDOM_117_pos2 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250825YNorth_test_table_Lab117Pos2.txt",header=0,sep=" ")
df_mDOM_117_pos3 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250825YNorth_test_table_Lab117Pos3.txt",header=0,sep=" ")
df_mDOM_117_floor = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250825YNorth_test_table_Lab117Floor.txt",header=0,sep=" ")

#measurement taken on Sep 8 2025
df_mDOM_rooftop_2 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250908XNorth_rooftop_near_telescope.txt",header=0,sep=" ")
df_mDOM_rooftop_away_2 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250908XNorth_rooftop_away_from_telescope.txt",header=0,sep=" ")
df_mDOM_117_2 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250908XNorth_Lab117.txt",header=0,sep=" ")

#measurement at mDOM at Utah on Oct 09 2025, PMT 8/16 facing North (the water tank), rotating anticlockwise when looking from top
df_mDOM_Utah = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/utah/mDOM_North_UtahScan_PMT8_16North.txt",header=0,sep=" ")
# df_mDOM_Utah = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/utah/mDOM_North_UtahScan.txt",header=0,sep=" ")


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

def get_mean_B_rolling_muT(df,angle):
    '''
    returns field in muT and degrees
    '''
    bx,by,bz = df[["bx","by","bz"]].values[angle:angle+10].T
    # print(f"bx {bx} by {by} bz {bz}")
    angle_list = []
    r,theta,phi = to_spherical_list(bx,by,bz)
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
    return np.mean(bx),np.std(bx),np.mean(by),np.std(by),np.mean(bz),np.std(bz),\
        np.mean(r),np.std(r),np.rad2deg(np.mean(theta)),np.rad2deg(np.std(theta)),np.rad2deg(np.mean(phi)),np.rad2deg(np.std(phi))

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
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
    ncolor = 0
    for df,df_label in zip(df_list,df_labels):
        angles_list = []
        r_list = []
        r_std_list = []
        theta_list = []
        theta_std_list = []
        phi_list = []
        phi_std_list = []
        bx_list = []
        bx_std_list = []
        by_list = []
        by_std_list = []
        bz_list = []
        bz_std_list = []
        for iangle in angles:
            bx,bx_std,by,by_std,bz,bz_std,r,r_std,theta,theta_std,phi,phi_std = get_mean_B_rolling(df,iangle)
            angles_list.append(iangle)
            bx_list.append(bx)
            bx_std_list.append(bx_std)
            by_list.append(by)
            by_std_list.append(by_std)
            bz_list.append(bz)
            bz_std_list.append(bz_std)
            r_list.append(r)
            r_std_list.append(r_std)
            theta_list.append(theta)
            theta_std_list.append(theta_std)
            phi_list.append(phi)
            phi_std_list.append(phi_std)
        bx_list_new = []
        by_list_new = []
        angles_list_new = []
        for ibx,iby,iangle in zip(bx_list,by_list,angles_list):
            if not np.isnan(ibx) and not np.isnan(iby):
                bx_list_new.append(ibx)
                by_list_new.append(iby)
                angles_list_new.append(iangle)
        bx_list = bx_list_new
        by_list = by_list_new
        angles_list = angles_list_new
        ax.plot(bx_list,by_list,"-o",c=colorsCustom[ncolor],label=f"{df_label}",alpha=1)
        bx_calibrated, by_calibrated = corrected_ellipse(bx_list, by_list)
        print(f"means bx {np.mean(bx_list)} {np.mean(bx_calibrated)}")
        print(f"means by {np.mean(by_list)} {np.mean(by_calibrated)}")
        ax.plot(bx_calibrated, by_calibrated, "--o",c=colorsCustom[ncolor+2], label=f"{df_label} (calibrated)", alpha=1)
        ncolor+=1
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xlim(0,360)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    ax.set_aspect('equal')
    ax.legend(loc="lower left",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" $B_x$ [$\mu$T]", fontsize=20)
    ax.set_ylabel(r" $B_y$ [$\mu$T]", fontsize=20)
    # ax1.set_ylim(35,130)
    # ax.set_ylim(-55,55)
    # ax.set_xlim(-55,55)
    # ax3.set_yticks(np.linspace(0,360,5))
    # ax.set_ylim(0,180)
    # ax2.set_ylim(-1.5,1.5)
    # ax3.set_ylim(8,11)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/orientation_with_B_rooftop_{MB}xy.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_B_rooftop_{MB}xy.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_orientation_360([df_mDOM_rooftop,df_mDOM_rooftop_pos2,df_mDOM_rooftop_away],["rooftop_near_telescope","rooftop_near_telescope_pos2","rooftop_away_telescope",],MB="mDOM")
# plot_orientation_360([df_mDOM_317,df_mMB_317],["mDOM MB 317", "mMB 317"],MB="mDOM_room")
# plot_xy_360([df_mDOM_117],["mDOM MB 117"],MB="mDOM_mb_117")
# plot_xy_360([df_mDOM_117,df_mDOM_117_repeat,df_mDOM_117_pos2,df_mDOM_117_pos3,df_mDOM_117_floor],["mDOM MB 117","mDOM MB 117 Repeat","mDOM MB 117 Pos2","mDOM MB 117 Pos3","mDOM MB 117 Floor"],MB="mDOM_lab")
# plot_xy_360([df_mDOM_117],["mDOM MB 117"],MB="mDOM_mb_117")
# plot_xy_360([df_mDOM_rooftop,df_mDOM_rooftop_pos2,df_mDOM_rooftop_away,df_mMB_rooftop],["rooftop_near_telescope","rooftop_near_telescope_pos2","rooftop_away_telescope","rooftop_near_telescope_mMB"],MB="mDOM_rooftop")
# plot_xy_360([df_mDOM_rooftop_2,df_mDOM_rooftop_away_2],["rooftop_near_telescope","rooftop_away_telescope"],MB="mDOM_rooftop_sep8")
# plot_xy_360([df_mDOM_rooftop_away],["rooftop_away_telescope",],MB="mDOM_rooftop")
# plot_xy_360([df_mDOM_Utah],["utah",],MB="mDOM_utah")

plot_xy_360([df_mDOM_117_2],["mDOM MB 117"],MB="mDOM_mb_117")


def plot_xyz_360(df_list,df_labels,MB):
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0],projection='3d')
    # for df,dir in zip(df_list,dir_list):
    ncolor = 0
    for df,df_label in zip(df_list,df_labels):
        angles_list = []
        r_list = []
        r_std_list = []
        theta_list = []
        theta_std_list = []
        phi_list = []
        phi_std_list = []
        bx_list = []
        bx_std_list = []
        by_list = []
        by_std_list = []
        bz_list = []
        bz_std_list = []
        for iangle in angles:
            bx,bx_std,by,by_std,bz,bz_std,r,r_std,theta,theta_std,phi,phi_std = get_mean_B_rolling(df,iangle)
            angles_list.append(iangle)
            bx_list.append(bx)
            bx_std_list.append(bx_std)
            by_list.append(by)
            by_std_list.append(by_std)
            bz_list.append(bz)
            bz_std_list.append(bz_std)
            r_list.append(r)
            r_std_list.append(r_std)
            theta_list.append(theta)
            theta_std_list.append(theta_std)
            phi_list.append(phi)
            phi_std_list.append(phi_std)
        bx_list_new = []
        by_list_new = []
        angles_list_new = []
        for ibx,iby,iangle in zip(bx_list,by_list,angles_list):
            if not np.isnan(ibx) and not np.isnan(iby):
                bx_list_new.append(ibx)
                by_list_new.append(iby)
                angles_list_new.append(iangle)
        bx_list = bx_list_new
        by_list = by_list_new
        angles_list = angles_list_new
        ax.plot3D(bx_list,by_list,bz_list,"-o",c=colorsCustom[ncolor],label=f"{df_label}",alpha=1)
        # bx_calibrated, by_calibrated = corrected_ellipse(bx_list, by_list)
        # ax.scatter3D(bx_calibrated, by_calibrated,bz_calibrated, "--o",c=colorsCustom[ncolor+2], label=f"{df_label} (calibrated)", alpha=1)
        ncolor+=1
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xlim(0,360)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    # ax.set_aspect('equal')
    ax.legend(loc="lower left",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" $B_x$ [$\mu$T]", fontsize=20)
    ax.set_ylabel(r" $B_y$ [$\mu$T]", fontsize=20)
    ax.set_zlabel(r" $B_z$ [$\mu$T]", fontsize=20)
    # ax1.set_ylim(35,130)
    # ax.set_ylim(-55,55)
    # ax.set_zlim(-60,0)
    # ax3.set_yticks(np.linspace(0,360,5))
    # ax.set_ylim(0,180)
    # ax2.set_ylim(-1.5,1.5)
    # ax3.set_ylim(8,11)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/orientation_with_B_rooftop_{MB}xyz.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_B_rooftop_{MB}xyz.pdf",transparent=False,bbox_inches='tight')
    # plt.close()
    plt.show()


# plot_xyz_360([df_mDOM_Utah],["utah",],MB="mDOM_utah")
plot_xyz_360([df_mDOM_117_2],["mDOM_117_2"],MB="mDOM_117")


mdom = "mDOM_MB_117"








def plot_corrected_heading_360(df_list,df_labels,MB,roll):
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
    ncolor = 0
    for df,df_label in zip(df_list,df_labels):
        angles_list = []
        r_list = []
        r_std_list = []
        theta_list = []
        theta_std_list = []
        phi_list = []
        phi_std_list = []
        bx_list = []
        bx_std_list = []
        by_list = []
        by_std_list = []
        bz_list = []
        bz_std_list = []
        for iangle in angles:
            bx,bx_std,by,by_std,bz,bz_std,r,r_std,theta,theta_std,phi,phi_std = get_mean_B_rolling(df,iangle)
            angles_list.append(iangle)
            bx_list.append(bx)
            bx_std_list.append(bx_std)
            by_list.append(by)
            by_std_list.append(by_std)
            bz_list.append(bz)
            bz_std_list.append(bz_std)
            r_list.append(r)
            r_std_list.append(r_std)
            theta_list.append(theta)
            theta_std_list.append(theta_std)
            phi_list.append(phi)
            phi_std_list.append(phi_std)
        # ax.plot(bx_list,by_list,"-o",label=f"{df_label}",alpha=0.5)
        bx_list_new = []
        by_list_new = []
        angles_list_new = []
        for ibx,iby,iangle in zip(bx_list,by_list,angles_list):
            if not np.isnan(ibx) and not np.isnan(iby):
                bx_list_new.append(ibx)
                by_list_new.append(iby)
                angles_list_new.append(iangle)
        bx_list = bx_list_new
        by_list = by_list_new
        angles_list = angles_list_new

        bx_calibrated, by_calibrated = corrected_ellipse(bx_list, by_list)
        headings_original = [np.rad2deg(np.arctan2(by, bx)) for bx, by in zip(bx_list, by_list)]
        print(f"headings original: {headings_original}")
        headings_original = [i+360 if i<-1.0 else i for i in headings_original]
        print(f"headings original: {headings_original}")
        headings_corrected = [np.rad2deg(np.arctan2(byc, bxc)) for bxc,byc in zip(bx_calibrated, by_calibrated)]
        print(f"headings corrected: {headings_corrected}")
        headings_corrected = [i+360 if i<-1.0 else i for i in headings_corrected]
        # headings_corrected = [i+360 if 0<i<42 else i for i in headings_corrected]
        headings_corrected = [i+360 if 0<i<34.0 else i for i in headings_corrected]
        # print(f"headings original: {headings_original}")
        print(f"headings corrected: {headings_corrected}")
        print(f"means bx {np.mean(bx_list)} {np.mean(bx_calibrated)}")
        print(f"means by {np.mean(by_list)} {np.mean(by_calibrated)}")
        # ax.plot(bx_calibrated, by_calibrated, "-o", label=f"{df_label} (calibrated)", alpha=0.5)
        # roll = True
        # if roll:
        #     angles_list = np.roll(angles_list, shift=11)
        #     headings_corrected = np.roll(headings_corrected, shift=11)
        #     headings_original = np.roll(headings_original, shift=11)
        #     if df_label == "mDOM MB 117":
        #         angles_list = [- (360 - i) if i >= 270 else i for i in angles_list]
        #     elif df_label in ["rooftop_away_telescope"]:
        #         angles_list = [- (360 - i) if i >= 290 else i for i in angles_list]
        #     elif df_label in ["rooftop_near_telescope_pos2"]:
        #         angles_list = [- (360 - i) if i >= 210 else i for i in angles_list]
        #     elif df_label in ["rooftop_near_telescope"]:
        #         angles_list = [- (360 - i) if i >= 220 else i for i in angles_list]
        #     elif df_label in ["rooftop_near_telescope_run2"]:
        #         angles_list = [- (360 - i) if i >= 330 else i for i in angles_list]
        #     elif df_label in ["rooftop_away_telescope_run2"]:
        #         angles_list = [- (360 - i) if i >= 330 else i for i in angles_list]
        #     else:
        #         angles_list = [- (360 - i) if i >= 260 else i for i in angles_list]
        # else:
        #     ax.plot(angles_list, headings_original, "-o", c=colorsCustom[ncolor], label=f"{df_label} (original)", alpha=1)
        ax.plot(angles_list, headings_original, "-o", c=colorsCustom[ncolor], label=f"{df_label}", alpha=1)
        ax.plot(angles_list, headings_corrected, "--o", c=colorsCustom[ncolor+2], label=f"{df_label} (calibrated)", alpha=1)
        ncolor += 1
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xlim(0,360)
    ###############################################
    ##############################################
    # ax.set_xticks(np.linspace(0,360,9)-roll)
    # ax.set_yticks(np.linspace(0,360,9))
    ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    ax.set_yticks(np.linspace(0,405,10))
    # ax.set_yticks(np.linspace(-180,180,9))
    #############################################
    #############################################
    ax.set_aspect('equal')
    ax.legend(ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" mDOM rotation [$^{\circ}$]", fontsize=20)
    ax.set_ylabel(r" $\phi$ [$^{\circ}$]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_B_rooftop_{MB}xy_heading_roll{roll}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_B_rooftop_{MB}xy_heading_roll{roll}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_orientation_360([df_mDOM_rooftop,df_mDOM_rooftop_pos2,df_mDOM_rooftop_away],["rooftop_near_telescope","rooftop_near_telescope_pos2","rooftop_away_telescope",],MB="mDOM")
# plot_orientation_360([df_mDOM_317,df_mMB_317],["mDOM MB 317", "mMB 317"],MB="mDOM_room")
# plot_xy_360([df_mDOM_117],["mDOM MB 117"],MB="mDOM_mb_117")
# plot_corrected_heading_360([df_mDOM_117_pos2],["mDOM MB 117 Pos2"],MB="mDOM_lab_117_pos2")
# plot_corrected_heading_360([df_mDOM_117,df_mDOM_117_repeat,df_mDOM_117_pos2,df_mDOM_117_pos3,df_mDOM_117_floor],["mDOM MB 117","mDOM MB 117 Repeat","mDOM MB 117 Pos2","mDOM MB 117 Pos3","mDOM MB 117 Floor"],MB="mDOM_lab")
# plot_corrected_heading_360([df_mDOM_117,df_mDOM_117_pos2,df_mDOM_117_pos3],["mDOM MB 117","mDOM MB 117 Pos2","mDOM MB 117 Pos3"],MB="mDOM_lab")
# plot_corrected_heading_360([df_mDOM_117_floor],["mDOM MB 117 Floor"],MB="mDOM_lab")
# plot_corrected_heading_360([df_mDOM_117,df_mDOM_117_repeat,df_mDOM_117_floor],["mDOM MB 117","mDOM MB 117 Repeat","mDOM MB 117 Floor"],MB="mDOM_lab")
# plot_xy_360([df_mDOM_117],["mDOM MB 117"],MB="mDOM_mb_117")
# plot_corrected_heading_360([df_mDOM_rooftop,df_mDOM_rooftop_pos2,df_mDOM_rooftop_away],["rooftop_near_telescope","rooftop_near_telescope_pos2","rooftop_away_telescope",],MB="mDOM_rooftop")
###################################################################################################################################################################
#############################latest measurement####################################################################################################################
# plot_corrected_heading_360([df_mDOM_rooftop_2,df_mDOM_rooftop_away_2],["rooftop_near_telescope_run2","rooftop_away_telescope_run2"],MB="mDOM_rooftop_sep8",roll=True)
# plot_corrected_heading_360([df_mDOM_rooftop_2,df_mDOM_rooftop_away_2],["rooftop_near_telescope_run2","rooftop_away_telescope_run2"],MB="mDOM_rooftop_sep8",roll=False)
# plot_xy_360([df_mDOM_117_2],["mDOM_MB_117"],MB="mDOM_117_sep8")
# plot_corrected_heading_360([df_mDOM_117_2],["mDOM_MB_117"],MB="mDOM_117_sep8",roll=45)    
plot_corrected_heading_360([df_mDOM_117_2],["mDOM_MB_117"],MB="mDOM_117_sep8",roll=0)
# plot_corrected_heading_360([df_mDOM_Utah],["mDOM_utah"],MB="mDOM_utah",roll=0)


def get_mean_B_list(df_list,df_labels):
    for df,df_label in zip(df_list,df_labels):
        angles_list = []
        r_list = []
        r_std_list = []
        theta_list = []
        theta_std_list = []
        phi_list = []
        phi_std_list = []
        bx_list = []
        bx_std_list = []
        by_list = []
        by_std_list = []
        bz_list = []
        bz_std_list = []
        for iangle in angles:
            bx,bx_std,by,by_std,bz,bz_std,r,r_std,theta,theta_std,phi,phi_std = get_mean_B_rolling(df,iangle)
            angles_list.append(iangle)
            bx_list.append(bx)
            bx_std_list.append(bx_std)
            by_list.append(by)
            by_std_list.append(by_std)
            bz_list.append(bz)
            bz_std_list.append(bz_std)
            r_list.append(r)
            r_std_list.append(r_std)
            theta_list.append(theta)
            theta_std_list.append(theta_std)
            phi_list.append(phi)
            phi_std_list.append(phi_std)
        # ax.plot(bx_list,by_list,"-o",label=f"{df_label}",alpha=0.5)
        bx_list_new = []
        by_list_new = []
        bz_list_new = []
        angles_list_new = []
        for ibx,iby,ibz,iangle in zip(bx_list,by_list,bz_list,angles_list):
            if not np.isnan(ibx) and not np.isnan(iby) and not np.isnan(ibz):
                bx_list_new.append(ibx)
                by_list_new.append(iby)
                bz_list_new.append(ibz)
                angles_list_new.append(iangle)
        bx_list = bx_list_new
        by_list = by_list_new
        bz_list = bz_list_new
        angles_list = angles_list_new
    return bx_list, by_list, bz_list, angles_list

bx_mean_list, by_mean_list, bz_mean_list, angles_list = get_mean_B_list([df_mDOM_117_2],["mDOM_MB_117"])

#at tuscaloosa
B_horizontal = 23.0 #muT
Bvert = 42.3 #muT
B_total = 48.1 #muT



def plot_B_calibrated(mag_x, mag_y, mag_z,label=""):
    mag_x = np.array(mag_x)*10**6 #convert to microTesla
    mag_y = np.array(mag_y)*10**6 #convert to microTesla
    mag_z = np.array(mag_z)*10**6 #convert to microTesla
    mag_x_cal, mag_y_cal = corrected_ellipse(mag_x, mag_y)
    B = np.sqrt(mag_x**2 + mag_y**2 + mag_z**2)
    B_cal = np.sqrt(np.array(mag_x_cal)**2 + np.array(mag_y_cal)**2 + mag_z**2)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    ax.plot(angles_list,B, "-o", c="b", label=f"{'raw'}", alpha=1)
    ax.plot(angles_list,B_cal, "-o", c="r", label=f"{'cal'}", alpha=1)
    B_cal_noaa = np.sqrt(np.array(mag_x_cal)**2 + np.array(mag_y_cal)**2 + (np.ones_like(mag_z)*Bvert)**2)
    # ax.plot(angles_list,B_cal_noaa, "-o", c="g", label=f"{'cal_noaa'}", alpha=1)
    ax.hlines(B_total,0,360,ls="--",lw=2.5,label=f"B$_{{total}}$ ({B_total:.1f} ${{\u03bc}}$T)",alpha=1.0)
    ax.set_xticks(np.linspace(0,360,9))
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.set_ylim(0,90)
    # ax.set_aspect('equal')
    ax.legend(loc="lower left",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_ylabel(r" $B$ [$\mu$T]", fontsize=20)
    ax.set_xlabel(r" $\phi$ [$^{\circ}$]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}B_calibrated_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}B_calibrated_{label}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_B_calibrated(bx_mean_list, by_mean_list, bz_mean_list,label="coarse")


def plot_Bxy_calibrated(mag_x, mag_y, mag_z,label=""):
    mag_x = np.array(mag_x)*10**6 #convert to microTesla
    mag_y = np.array(mag_y)*10**6 #convert to microTesla
    mag_z = np.array(mag_z)*10**6 #convert to microTesla
    mag_x_cal, mag_y_cal = corrected_ellipse(mag_x, mag_y)
    B = np.sqrt(mag_x**2 + mag_y**2)
    B_cal = np.sqrt(np.array(mag_x_cal)**2 + np.array(mag_y_cal)**2)
    print(f"B-Bcal range: {np.min(abs(B-B_cal)):.2f} to {np.max(abs(B-B_cal)):.2f} microTesla")
    print(f"B-Bcal range: {[i for i in B-B_cal if i < 0]} microTesla")
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    ax.hlines(B_horizontal,0,360,ls="--",lw=2.5,label=f"B$_{{xy}}$ ({B_horizontal:.1f} ${{\u03bc}}$T)",alpha=1.0)
    ax.plot(angles_list,B, "-o", c="b", label=f"{'raw'}", alpha=1)
    ax.plot(angles_list,B_cal, "-o", c="r", label=f"{'cal'}", alpha=1)
    ax.set_xticks(np.linspace(0,360,9))
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.set_aspect('equal')
    ax.set_ylim(0,90)
    ax.legend(loc="upper right",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_ylabel(r" $B_{xy}$ [$\mu$T]", fontsize=20)
    ax.set_xlabel(r" $\phi$ [$^{\circ}$]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}BxBy_calibrated_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}BxBy_calibrated_{label}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_Bxy_calibrated(bx_mean_list, by_mean_list, bz_mean_list,label="coarse")




def plot_Bz(mag_z,label=""):
    mag_z = np.array(mag_z)*10**6 #convert to microTesla
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    ax.hlines(Bvert,0,360,ls="--",lw=2.5,label=f"B$_{{z}}$ ({Bvert:.1f} ${{\u03bc}}$T)",alpha=1.0)
    ax.plot(angles_list,abs(mag_z), "-o", c="b", label=f"{'raw'}", alpha=1)
    ax.set_xticks(np.linspace(0,360,9))
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    print(f"mean field offset Bz value: {Bvert-np.mean(mag_z):.2f} microTesla")
    print(f"range of Bz values: {np.min(Bvert-mag_z):.2f} to {np.max(Bvert-mag_z):.2f} microTesla")
    # ax.set_aspect('equal')
    ax.legend(loc="upper right",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_ylim(0,90)
    ax.set_ylabel(r" $B_{z}$ [$\mu$T]", fontsize=20)
    ax.set_xlabel(r" $\phi$ [$^{\circ}$]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}Bz_calibrated_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}Bz_calibrated_{label}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_Bz(bz_mean_list,label="coarse")
