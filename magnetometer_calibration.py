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





# print(df_mDOMMB_yNorth.head())

def get_mean_B(df):
    '''
    returns field in muT and degrees
    '''
    bx,by,bz = df[["bx","by","bz"]].values.T
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
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
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
        ax.plot(bx_list,by_list,"-o",label=f"{df_label}",alpha=0.5)
        bx_calibrated, by_calibrated = corrected_ellipse(bx_list, by_list)
        print(f"means bx {np.mean(bx_list)} {np.mean(bx_calibrated)}")
        print(f"means by {np.mean(by_list)} {np.mean(by_calibrated)}")
        ax.plot(bx_calibrated, by_calibrated, "-o", label=f"{df_label} (calibrated)", alpha=0.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xlim(0,360)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    ax.set_aspect('equal')
    ax.legend(ncols=2,fontsize=8)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" x [$\mu$T]", fontsize=16)
    ax.set_ylabel(r" y [$\mu$T]", fontsize=16)
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
plot_xy_360([df_mDOM_117,df_mDOM_117_repeat,df_mDOM_117_pos2,df_mDOM_117_pos3,df_mDOM_117_floor],["mDOM MB 117","mDOM MB 117 Repeat","mDOM MB 117 Pos2","mDOM MB 117 Pos3","mDOM MB 117 Floor"],MB="mDOM_lab")
# plot_xy_360([df_mDOM_117],["mDOM MB 117"],MB="mDOM_mb_117")
# plot_orientation_360([df_mDOM_rooftop,df_mDOM_rooftop_pos2,df_mDOM_rooftop_away,df_mMB_rooftop],["rooftop_near_telescope","rooftop_near_telescope_pos2","rooftop_away_telescope","rooftop_near_telescope_mMB"],MB="mDOM")



def plot_corrected_heading_360(df_list,df_labels,MB):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
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
        bx_calibrated, by_calibrated = corrected_ellipse(bx_list, by_list)
        headings_original = [np.rad2deg(np.arctan2(by, bx)) for bx, by in zip(bx_list, by_list)]
        headings_original = [i+360 if i<0 else i for i in headings_original]
        headings_corrected = [np.rad2deg(np.arctan2(byc, bxc)) for bxc,byc in zip(bx_calibrated, by_calibrated)]
        headings_corrected = [i+360 if i<0 else i for i in headings_corrected]
        print(f"means bx {np.mean(bx_list)} {np.mean(bx_calibrated)}")
        print(f"means by {np.mean(by_list)} {np.mean(by_calibrated)}")
        # ax.plot(bx_calibrated, by_calibrated, "-o", label=f"{df_label} (calibrated)", alpha=0.5)
        ax.plot(angles, headings_original, "-o", label=f"{df_label} (original)", alpha=0.5)
        ax.plot(angles, headings_corrected, "-o", label=f"{df_label} (calibrated)", alpha=0.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xlim(0,360)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    # ax.set_aspect('equal')
    ax.legend(ncols=2,fontsize=8)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" x [$\mu$T]", fontsize=16)
    ax.set_ylabel(r" y [$\mu$T]", fontsize=16)
    plt.savefig(plotFolder+f"/orientation_with_B_rooftop_{MB}xy_heading.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_B_rooftop_{MB}xy_heading.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_orientation_360([df_mDOM_rooftop,df_mDOM_rooftop_pos2,df_mDOM_rooftop_away],["rooftop_near_telescope","rooftop_near_telescope_pos2","rooftop_away_telescope",],MB="mDOM")
# plot_orientation_360([df_mDOM_317,df_mMB_317],["mDOM MB 317", "mMB 317"],MB="mDOM_room")
# plot_xy_360([df_mDOM_117],["mDOM MB 117"],MB="mDOM_mb_117")
# plot_corrected_heading_360([df_mDOM_117_pos2],["mDOM MB 117 Pos2"],MB="mDOM_lab_117_pos2")
# plot_corrected_heading_360([df_mDOM_117,df_mDOM_117_repeat,df_mDOM_117_pos2,df_mDOM_117_pos3,df_mDOM_117_floor],["mDOM MB 117","mDOM MB 117 Repeat","mDOM MB 117 Pos2","mDOM MB 117 Pos3","mDOM MB 117 Floor"],MB="mDOM_lab")
plot_corrected_heading_360([df_mDOM_117,df_mDOM_117_pos2,df_mDOM_117_pos3],["mDOM MB 117","mDOM MB 117 Pos2","mDOM MB 117 Pos3"],MB="mDOM_lab")
# plot_corrected_heading_360([df_mDOM_117_floor],["mDOM MB 117 Floor"],MB="mDOM_lab")
# plot_corrected_heading_360([df_mDOM_117,df_mDOM_117_repeat,df_mDOM_117_floor],["mDOM MB 117","mDOM MB 117 Repeat","mDOM MB 117 Floor"],MB="mDOM_lab")
# plot_xy_360([df_mDOM_117],["mDOM MB 117"],MB="mDOM_mb_117")
# plot_orientation_360([df_mDOM_rooftop,df_mDOM_rooftop_pos2,df_mDOM_rooftop_away,df_mMB_rooftop],["rooftop_near_telescope","rooftop_near_telescope_pos2","rooftop_away_telescope","rooftop_near_telescope_mMB"],MB="mDOM")
