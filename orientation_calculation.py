#!/usr/bin/env python3

import os
import glob

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from utils import to_spherical, tilt_angle, get_means, get_sub_df, to_spherical_list, to_spherical_list_rotated, meas_from_df

from pathlib import Path
home = str(Path.home())

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = home+"/research_ua/icecube/upgrade/upgrade_comissioning/plots/"
# print(degg_list)

colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']
B_gallallee = 48.1428 #muT


drive = "/Users/epaudel/Library/CloudStorage/OneDrive-TheUniversityofAlabama/"
accelerometer_readings = drive+"research_notes/references/upgrade_comissioning/mDOM_orientations/accelerometer_measurement.xlsx"

#new measurement taken on May 29 2025 for full 360 rotation and acceclerometer mMB and mDOMMB
mDOM_tilt_data = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/mDOMMB_accelerometer_magnetometerTiltMay29.txt"
mMB_tilt_data = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/mMB_accelerometer_magnetometerTiltMay29.txt"


import pandas as pd

df_mDOMMB = pd.read_csv(mDOM_tilt_data,header=0,sep=" ")
df_mMB = pd.read_csv(mMB_tilt_data,header=0,sep=" ")

measured_angles = [f"{i}{j}" for i in range(0,361,10) for j in ["","A","B","C","D"]]
print(measured_angles)



accelerometer_cart_dict_mDOMMB,accelerometer_sphe_dict_mDOMMB,magnetometer_cart_dict_mDOMMB,magnetometer_sphe_dict_mDOMMB = meas_from_df(df_mDOMMB,measured_angles,sensor_rotation=False)
# accelerometer_cart_dict_mMB,accelerometer_sphe_dict_mMB,magnetometer_cart_dict_mMB,magnetometer_sphe_dict_mMB = meas_from_df(df_mMB,sensor_rotation=True) #setup 1 May 16
accelerometer_cart_dict_mMB,accelerometer_sphe_dict_mMB,magnetometer_cart_dict_mMB,magnetometer_sphe_dict_mMB = meas_from_df(df_mMB,measured_angles,sensor_rotation=False) #setup 2 May 29

plane = [f"{n}" for n in range(0,361,10)]
A_tilts = [ix for ix in measured_angles if "A" in ix]
B_tilts = [ix for ix in measured_angles if "B" in ix]
C_tilts = [ix for ix in measured_angles if "C" in ix]
D_tilts = [ix for ix in measured_angles if "D" in ix]

plane_int = [int(i) for i in range(0,361,10)]

# rotation matrices
def rotation_matrix_x(theta):
    return np.array([[1, 0, 0],
                     [0, np.cos(theta), -np.sin(theta)],
                     [0, np.sin(theta), np.cos(theta)]])

def rotation_matrix_y(theta):
    return np.array([[np.cos(theta), 0, np.sin(theta)],
                     [0, 1, 0],
                     [-np.sin(theta), 0, np.cos(theta)]])

def rotation_matrix_z(theta):
    return np.array([[np.cos(theta), -np.sin(theta), 0],
                     [np.sin(theta), np.cos(theta), 0],
                     [0, 0, 1]])

def calculate_heading(gx,gy,gz,bx,by,bz):
    # normalize accelerometer data
    g = np.sqrt(gx**2 + gy**2 + gz**2)
    gx /= g
    gy /= g
    gz /= g

    #theta and phi angles from accelerometer
    theta = np.arctan2(np.sqrt(gx**2 + gy**2), gz)
    phi = np.arctan2(gy, gx)

    #correct magnetometer readings for tilt
    # bx_rotated,by_rotated,bz_rotated = np.dot(rotation_matrix_x(theta), np.array([bx, by, bz]))
    bx_rotated,by_rotated,bz_rotated = np.dot(rotation_matrix_z(phi),np.dot(rotation_matrix_x(theta), np.array([bx, by, bz])))
    heading = np.arctan2(by_rotated, bx_rotated)
    b_tilt = np.arctan2(np.sqrt(bx_rotated**2 + by_rotated**2), bz_rotated)
    if heading < 0:
        heading += 2 * np.pi

    return heading,b_tilt


def plot_sensor_heading(accelerometer_cart_dict,accelerometer_sphe_dict,magnetometer_cart_dict,magnetometer_sphe_dict,MB):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0])
    # for ilabel,isetting in zip(["no tilt"],[plane]):
    for ilabel,isetting in zip(["no tilt","tab A","tab B","tab C","tab D"],[plane, A_tilts, B_tilts,C_tilts,D_tilts]):
        angles = []
        b_list = []
        b_err_list = []
        theta_list = []
        theta_err_list = []
        phi_list = []
        phi_err_list = []
        heading_list = []
        b_tilt_list = []
        for iangle in isetting:
            bx,bx_err,by,by_err,bz,bz_err = magnetometer_cart_dict[iangle]
            gx,gx_err,gy,gy_err,gz,gz_err = accelerometer_cart_dict[iangle]
            b,b_err,theta,theta_err,phi,phi_err = magnetometer_sphe_dict[iangle]
            theta_list.append(np.rad2deg(theta))
            phi_list.append(np.rad2deg(phi))
            heading,b_tilt = calculate_heading(gx,gy,gz,bx,by,bz)
            heading_list.append(np.rad2deg(heading))
            b_tilt_list.append(np.rad2deg(b_tilt))
        ax.plot(plane_int,phi_list,"o-",label=f"{ilabel}_phi")
        ax.plot(plane_int,heading_list,"<-",label=f"{ilabel}_heading")
        # ax.plot(plane_int,theta_list,"^-",label=f"{ilabel}_theta")
        # ax.plot(plane_int,b_tilt_list,"s-",label=f"{ilabel}_b_tilt")

        ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
        ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)        
        ax.grid(True,alpha=0.6)
        # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
        # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
        # ax.set_xlim(0,360)
        ax.set_xticks(np.linspace(0,360,9))
        # ax.set_yticks(np.linspace(0,360,9))
        # ax.set_aspect('equal')
        ax.legend(ncols=5)    
    ax.set_ylabel(r" $\phi$ [$^{\circ}$]", fontsize=22)
    # ax2.set_ylim(-1.5,1.5)
    # ax3.set_ylim(8,11)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/sensor_orientationg{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/sensor_orientationg{MB}.pdf",transparent=False,bbox_inches='tight')
    plt.close()    

plot_sensor_heading(accelerometer_cart_dict_mDOMMB,accelerometer_sphe_dict_mDOMMB,magnetometer_cart_dict_mDOMMB,magnetometer_sphe_dict_mDOMMB,MB="mDOM")