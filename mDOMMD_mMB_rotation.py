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
 
# mDOM_tilt_data = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/mDOMMB_accelerometer_magnetometerTiltMay16.txt"
# mMB_tilt_data = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/mMB_accelerometer_magnetometerTiltMay16.txt"
#data format {gx} {gy} {gz} {bx} {by} {bz} for 0 to 90 degree angles and ABCD tilts

#new measurement taken on May 29 2025 for full 360 rotation and acceclerometer mMB and mDOMMB
mDOM_tilt_data = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/mDOMMB_accelerometer_magnetometerTiltMay29.txt"
mMB_tilt_data = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/mMB_accelerometer_magnetometerTiltMay29.txt"



import pandas as pd

df_mDOMMB = pd.read_csv(mDOM_tilt_data,header=0,sep=" ")
df_mMB = pd.read_csv(mMB_tilt_data,header=0,sep=" ")
# print(df.head)
# print(df.columns.values)
# print(df[df["angle"]== "0"])
# measured_angles = list(set(df_mDOMMB["angle"].values))
# print(f"measured angles {measured_angles}")
# measured_angles = sorted(measured_angles)
measured_angles = [f"{i}{j}" for i in range(0,361,10) for j in ["","A","B","C","D"]]
print(measured_angles)
# print(f"measured angles {measured_angles}")

accelerometer_cart_dict_mDOMMB,accelerometer_sphe_dict_mDOMMB,magnetometer_cart_dict_mDOMMB,magnetometer_sphe_dict_mDOMMB = meas_from_df(df_mDOMMB,measured_angles,sensor_rotation=False)
# accelerometer_cart_dict_mMB,accelerometer_sphe_dict_mMB,magnetometer_cart_dict_mMB,magnetometer_sphe_dict_mMB = meas_from_df(df_mMB,sensor_rotation=True) #setup 1 May 16
accelerometer_cart_dict_mMB,accelerometer_sphe_dict_mMB,magnetometer_cart_dict_mMB,magnetometer_sphe_dict_mMB = meas_from_df(df_mMB,measured_angles,sensor_rotation=False) #setup 2 May 29

# plane = [f"{n}" for n in range(0,91,10)]
plane = [f"{n}" for n in range(0,361,10)]
A_tilts = [ix for ix in measured_angles if "A" in ix]
B_tilts = [ix for ix in measured_angles if "B" in ix]
C_tilts = [ix for ix in measured_angles if "C" in ix]
D_tilts = [ix for ix in measured_angles if "D" in ix]

plane_int = [int(i) for i in range(0,361,10)]

def plot_accelerometer_xyz(accelerometer_cart_dict,MB):
    fig = plt.figure(figsize=(8,5*3))
    gs = gridspec.GridSpec(nrows=3,ncols=1, figure=fig)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[2,0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    # print(len(rotations_list),len(azimuth),azimuth)
    for ilabel,isetting in zip(["no tilt","tab A","tab B","tab C","tab D"],[plane, A_tilts, B_tilts,C_tilts,D_tilts]):
        angles = []
        gx_list = []
        gx_err_list = []
        gy_list = []
        gy_err_list = []
        gz_list = []
        gz_err_list = []
        for iangle in isetting:
            print("############################################################################")
            print(iangle)
            # print(accelerometer_cart_dict)
            print(accelerometer_cart_dict[iangle])
            gx,gx_err,gy,gy_err,gz,gz_err = accelerometer_cart_dict[iangle]
            gx_list.append(gx)
            gx_err_list.append(gx_err)
            gy_list.append(gy)
            gy_err_list.append(gy_err)
            gz_list.append(gz)
            gz_err_list.append(gz_err)
        ax1.errorbar(plane_int,gx_list,yerr=gx_err_list,fmt="o",label=f"{ilabel}")
        ax2.errorbar(plane_int,gy_list,yerr=gy_err_list,fmt="o",label=f"{ilabel}")
        ax3.errorbar(plane_int,gz_list,yerr=gz_err_list,fmt="o",label=f"{ilabel}")


        # nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        # ax.errorbar(rotations_list,gxy_means,yerr=gxy_stds,fmt="o",label=f"{isetting}")
    # ax.errorbar(rotations,gxy_means,yerr=gxy_stds,fmt="o",label="azimuth")
    # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
    for ax in [ax1,ax2,ax3]:
        ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
        ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
        
        ax.grid(True,alpha=0.6)
        ax.set_xticks(np.linspace(0,360,9))
        # ax.set_yticks(np.linspace(0,360,9))
        # ax.set_aspect('equal')
        ax.legend(ncols=5)
    ax1.set_ylabel(r" g$_x$ [ms$^{-2}$]", fontsize=22)
    ax2.set_ylabel(r" g$_y$ [ms$^{-2}$]", fontsize=22)
    ax3.set_ylabel(r" g$_z$ [ms$^{-2}$]", fontsize=22)
    # ax1.set_ylim(-1.5,1.5)
    # ax2.set_ylim(-1.5,1.5)
    # ax3.set_ylim(8,11)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/magnetometerAccelerometerTilt_gxgygz{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometerAccelerometerTilt_gxgygz{MB}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_accelerometer_xyz(accelerometer_cart_dict_mDOMMB,MB="mDOMMB")
plot_accelerometer_xyz(accelerometer_cart_dict_mMB,MB="mMB")


def plot_accelerometer_rθφ(accelerometer_sphe_dict,MB):
    fig = plt.figure(figsize=(8,5*3))
    gs = gridspec.GridSpec(nrows=3,ncols=1, figure=fig)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[2,0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    # print(len(rotations_list),len(azimuth),azimuth)
    for ilabel,isetting in zip(["no tilt","tab A","tab B","tab C","tab D"],[plane, A_tilts, B_tilts,C_tilts,D_tilts]):
        angles = []
        g_list = []
        g_err_list = []
        theta_list = []
        theta_err_list = []
        phi_list = []
        phi_err_list = []
        for iangle in isetting:
            g,g_err,theta,theta_err,phi,phi_err = accelerometer_sphe_dict[iangle]
            g_list.append(g)
            g_err_list.append(g_err)
            theta_list.append(np.rad2deg(theta))
            theta_err_list.append(np.rad2deg(theta_err))
            phi_list.append(np.rad2deg(phi))
            phi_err_list.append(np.rad2deg(phi_err))
        ax1.errorbar(plane_int,g_list,yerr=g_err_list,fmt="o",label=f"{ilabel}")
        ax2.errorbar(plane_int,theta_list,yerr=theta_err_list,fmt="o",label=f"{ilabel}")
        ax3.errorbar(plane_int,phi_list,yerr=phi_err_list,fmt="o",label=f"{ilabel}")


        # nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        # ax.errorbar(rotations_list,gxy_means,yerr=gxy_stds,fmt="o",label=f"{isetting}")
    # ax.errorbar(rotations,gxy_means,yerr=gxy_stds,fmt="o",label="azimuth")
    # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
    for ax in [ax1,ax2,ax3]:
        ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
        ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)        
        ax.grid(True,alpha=0.6)
        ax.set_xticks(np.linspace(0,360,9))
        # ax.set_yticks(np.linspace(0,360,9))
        # ax.set_aspect('equal')
        ax.legend(ncols=5)
    ax1.set_ylabel(r" g [ms$^{-2}$]", fontsize=22)
    ax1.set_ylim(9.7,9.915)
    ax2.set_ylabel(r" $\theta$ [$^{\circ}$]", fontsize=22)
    ax3.set_ylabel(r" $\phi$ [$^{\circ}$]", fontsize=22)
    ax3.set_yticks(np.linspace(0,360,5))
    # ax1.set_ylim(-1.5,1.5)
    # ax2.set_ylim(-1.5,1.5)
    # ax3.set_ylim(8,11)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/magnetometerAccelerometerTilt_gthetaphi{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometerAccelerometerTilt_gthetaphi{MB}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_accelerometer_rθφ(accelerometer_sphe_dict_mDOMMB,MB="mDOMMB")
plot_accelerometer_rθφ(accelerometer_sphe_dict_mMB,MB="mMB")

def plot_magnetometer_xyz(magnetometer_cart_dict,MB):
    fig = plt.figure(figsize=(8,5*3))
    gs = gridspec.GridSpec(nrows=3,ncols=1, figure=fig)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[2,0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    # print(len(rotations_list),len(azimuth),azimuth)
    for ilabel,isetting in zip(["no tilt","tab A","tab B","tab C","tab D"],[plane, A_tilts, B_tilts,C_tilts,D_tilts]):
        angles = []
        bx_list = []
        bx_err_list = []
        by_list = []
        by_err_list = []
        bz_list = []
        bz_err_list = []
        for iangle in isetting:
            bx,bx_err,by,by_err,bz,bz_err = magnetometer_cart_dict[iangle]
            bx_list.append(bx*10**(6)) #scaling T --> muT
            bx_err_list.append(bx_err*10**(6))
            by_list.append(by*10**(6))
            by_err_list.append(by_err*10**(6))
            bz_list.append(bz*10**(6))
            bz_err_list.append(bz_err*10**(6))
        ax1.errorbar(plane_int,bx_list,yerr=bx_err_list,fmt="o",label=f"{ilabel}")
        ax2.errorbar(plane_int,by_list,yerr=by_err_list,fmt="o",label=f"{ilabel}")
        ax3.errorbar(plane_int,bz_list,yerr=bz_err_list,fmt="o",label=f"{ilabel}")


        # nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        # ax.errorbar(rotations_list,gxy_means,yerr=gxy_stds,fmt="o",label=f"{isetting}")
    # ax.errorbar(rotations,gxy_means,yerr=gxy_stds,fmt="o",label="azimuth")
    # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
    for ax in [ax1,ax2,ax3]:
        ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
        ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
        
        ax.grid(True,alpha=0.6)
        ax.set_xticks(np.linspace(0,360,9))
        # ax.set_yticks(np.linspace(0,360,9))
        # ax.set_aspect('equal')
        ax.legend(ncols=5)
    ax1.set_ylabel(r" B$_x$ [$\mu$T]", fontsize=22)
    ax2.set_ylabel(r" B$_y$ [$\mu$T]", fontsize=22)
    ax3.set_ylabel(r" B$_z$ [$\mu$T]", fontsize=22)
    # ax1.set_ylim(-1.5,1.5)
    # ax2.set_ylim(-1.5,1.5)
    # ax3.set_ylim(8,11)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/magnetometerAccelerometerTilt_bxbybz{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometerAccelerometerTilt_bxbybz{MB}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_magnetometer_xyz(magnetometer_cart_dict_mDOMMB,MB="mDOMMB")
plot_magnetometer_xyz(magnetometer_cart_dict_mMB,MB="mMB")


def plot_magnetometer_rθφ(magnetometer_sphe_dict,MB):
    fig = plt.figure(figsize=(8,5*3))
    gs = gridspec.GridSpec(nrows=3,ncols=1, figure=fig)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[2,0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    # print(len(rotations_list),len(azimuth),azimuth)
    for ilabel,isetting in zip(["no tilt","tab A","tab B","tab C","tab D"][:1],[plane, A_tilts, B_tilts,C_tilts,D_tilts][:1]):
        angles = []
        b_list = []
        b_err_list = []
        theta_list = []
        theta_err_list = []
        phi_list = []
        phi_err_list = []
        for iangle in isetting:
            b,b_err,theta,theta_err,phi,phi_err = magnetometer_sphe_dict[iangle]
            b_list.append(b*10**(6))
            b_err_list.append(b_err*10**(6))
            theta_list.append(np.rad2deg(theta))
            theta_err_list.append(np.rad2deg(theta_err))
            phi_list.append(np.rad2deg(phi))
            phi_err_list.append(np.rad2deg(phi_err))
        ax1.errorbar(plane_int,b_list,yerr=b_err_list,fmt="o",label=f"{ilabel}")
        
        ax2.errorbar(plane_int,theta_list,yerr=theta_err_list,fmt="o",label=f"{ilabel}")
        ax3.errorbar(plane_int,phi_list,yerr=phi_err_list,fmt="o",label=f"{ilabel}")


        # nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        # ax.errorbar(rotations_list,gxy_means,yerr=gxy_stds,fmt="o",label=f"{isetting}")
    # ax.errorbar(rotations,gxy_means,yerr=gxy_stds,fmt="o",label="azimuth")
    # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
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
        ax.set_xticks(np.linspace(0,360,9))
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
    plt.savefig(plotFolder+f"/magnetometerAccelerometerTilt_bthetaphi{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/magnetometerAccelerometerTilt_bthetaphi{MB}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_magnetometer_rθφ(magnetometer_sphe_dict_mDOMMB,MB="mDOMMB")
plot_magnetometer_rθφ(magnetometer_sphe_dict_mMB,MB="mMB")













# df = pd.read_excel(accelerometer_readings,sheet_name='accelerometer_05132025MiniMBDup',header=0,usecols="A:AK")#Aup arrow tab up by 44+-1 mm
# df_tilts = pd.read_excel(accelerometer_readings,header=0,usecols="R:V")
# rotations = df.columns.values
# tilts = df_tilts.columns.values
# print(df.columns.values)

# def get_column_names(name_string):
#     return [iname for iname in df.columns.values if name_string in iname]




# def get_means(df):
#     g_means = []
#     g_stds = []
#     gxy_means = []
#     gxy_stds = []
#     zenith_means = []
#     zenith_stds = []
#     azimuth_means = []
#     azimuth_stds = []
#     nx_means = []
#     ny_means = []
#     nz_means = []
#     nx_stds = []
#     ny_stds = []
#     nz_stds = []
#     for irot in df.columns.values[:]:
#         # print(f"{irot} rotation")
#         accelerometer_readings = df[[irot]].values    
#         # print("accelerometer reading",accelerometer_readings)
#         angles = []
#         zenith = []
#         azimuth = []
#         g = []
#         gxy = []
#         nx_list = []
#         ny_list = []
#         nz_list = []
#         for imeas in accelerometer_readings[:]:
#             # print(imeas[0])
#             # nx,ny,nz = imeas[0]
#             # print(imeas[0].split(","))
#             nxp,nyp,nzp = [float(jx.strip("[").strip("]")) for jx in imeas[0].split(",")]
#             nz = nxp
#             nx = nzp
#             ny = -nyp
#             # print(tilt_angle(nx,ny,nz)*180/np.pi)
#             angles.append(tilt_angle(nx,ny,nz)*180/np.pi)
#             zenith.append(to_spherical(nx,ny,nz)[1]*180/np.pi)
#             azimuth.append(to_spherical(nx,ny,nz)[2]*180/np.pi)
#             g.append(to_spherical(nx,ny,nz)[0]*10**6)
#             gxy.append(np.sqrt(nx**2+ny**2))
#             nx_list.append(nx)
#             ny_list.append(ny)
#             nz_list.append(nz)
#             # print(f"azimuth {azimuth}")
#         # print(f"Rotation {irot} mean {np.mean(angles):.2f} std {np.std(angles):.2f}")s
#         # print(f"Rotation {irot} mean zen {np.mean(zenith):.2f} std {np.std(zenith):.2f}"+
#         #       f" mean azi {np.mean(azimuth):.2f} std {np.std(azimuth):.2f}"+
#         #       f" mean B {np.mean(B):.2f} std {np.std(B):.2f}")
#         azimuth = [i+360 if i<0 else i for i in azimuth]
#         azimuth_means.append(np.mean(azimuth))
#         azimuth_stds.append(np.std(azimuth))
#         zenith_means.append(np.mean(zenith))
#         zenith_stds.append(np.std(azimuth))
#         g_means.append(np.mean(g))
#         g_stds.append(np.std(g))
#         gxy_means.append(np.mean(gxy))
#         gxy_stds.append(np.std(gxy))
#         nx_means.append(np.mean(nx_list))
#         ny_means.append(np.mean(ny_list))
#         nz_means.append(np.mean(nz_list))
#         nx_stds.append(np.std(nx_list))
#         ny_stds.append(np.std(ny_list))
#         nz_stds.append(np.std(nz_list))

#     # print(azimuth_means,g_means)
#     return nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds
# # nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(df_B_up)
# # print(df.iloc[:,:1].values)
# tab_labels = {"plane":"plane","A_up":"tab A","B_up":"tab B","C_up":"tab C","D_up":"tab D"}
# def plot_accelerometer_variation(df,rotation_list):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
#     # print(len(rotations_list),len(azimuth),azimuth)

#     _,_,_,_,_,_,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(df)
#     ax.errorbar(rotations,azimuth_means,yerr=azimuth_stds,fmt="o",label=f"")
#     # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
#     ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
#     ax.set_ylabel(r" azimuth angle of $\vec{g}$ [$^{\circ}$]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     ax.set_xticks(np.linspace(0,360,9))
#     ax.set_yticks(np.linspace(0,360,9))
#     # ax.set_aspect('equal')
#     ax.legend(ncols=4)
#     ax.set_xlim(0,360)
#     ax.set_ylim(0,360)
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/accelerometerxyz_fine.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/accelerometerxyz_fine.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# plot_accelerometer_variation(df,rotations)

# def plot_accelerometer_variation_zenith(df,rotation_list):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
#     # print(len(rotations_list),len(azimuth),azimuth)
#     _,_,_,_,_,_,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(df)
#     print(f"#######################zenith means {min(zenith_means)} {max(zenith_means)} {np.mean(zenith_means)}")
#     ax.errorbar(rotations,zenith_means,yerr=zenith_stds,fmt="o",label=f"")
#     # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
#     ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
#     ax.set_ylabel(r" zenith angle of $\vec{g}$ [$^{\circ}$]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     ax.set_xticks(np.linspace(0,360,9))
#     # ax.set_yticks(np.linspace(0,360,9))
#     # ax.set_aspect('equal')
#     ax.legend()
#     ax.set_xlim(0,360)
#     # ax.set_ylim(0,1.25)
#     # ax.set_ylim(0,12)
#     ax.set_ylim(0,3)
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/accelerometerxyz_fineZenith.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/accelerometerxyz_fineZenith.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# plot_accelerometer_variation_zenith(df,rotations)



def plot_gxy(gxy_means, Bx_stds):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    # print(len(rotations_list),len(azimuth),azimuth)
    for isetting in ["plane","A_up","B_up","C_up","D_up"]:
        this_df = df[get_column_names(isetting)]
        nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        ax.errorbar(rotations_list,gxy_means,yerr=gxy_stds,fmt="o",label=f"{isetting}")
    # ax.errorbar(rotations,gxy_means,yerr=gxy_stds,fmt="o",label="azimuth")
    # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
    ax.set_ylabel(r" gxy[$\mu$T$]", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    # ax.set_aspect('equal')
    ax.legend()
    # ax.set_xlim(0,360)
    # ax.set_ylim(0,360)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/accelerometergxy_fine.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/accelerometergxy_fine.pdf",transparent=False,bbox_inches='tight')
    plt.close()
# plot_gxy(gxy_means,gxy_stds)



def plot_accelerometer_variation_xyz(df,rotation_list):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    # print(len(rotations_list),len(azimuth),azimuth)
    for isetting in ["plane","A_up","B_up","C_up","D_up"][:1]:
        this_df = df[get_column_names(isetting)]
        nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        # ax.errorbar(rotations_list,nx_means,yerr=nx_stds,fmt="o",label=f"{isetting}"+"x")
        # ax.errorbar(rotations_list,ny_means,yerr=ny_stds,fmt="<",label=f"{isetting}"+"y")
        ax.errorbar(rotations_list,nz_means,yerr=nz_stds,fmt=">",label=f"g$_{{z}}$ {isetting}")
    # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
    # ax.set_ylabel(r" g_{x}", fontsize=22)
    # ax.set_ylabel(r" g_{y}", fontsize=22)
    ax.set_ylabel(r" g$_{z}$ [ms$^{-2}$]", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    # ax.set_aspect('equal')
    ax.legend()
    # ax.set_xlim(0,360)
    # ax.set_ylim(0,10)
    # ax.set_yscale("log")
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/accelerometerxyz_finenxnynz.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/accelerometerxyz_finenxnynz.pdf",transparent=False,bbox_inches='tight')
    plt.close()
# plot_accelerometer_variation_xyz(df,rotations_list)
axis_label = {"x"}
def plot_accelerometer_variation_xyz(df,rotation_list,which_axis):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    # print(len(rotations_list),len(azimuth),azimuth)
    for isetting in ["plane","A_up","B_up","C_up","D_up"][1:2]:
        this_df = df[get_column_names(isetting)]
        nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        # ax.errorbar(rotations_list,nx_means,yerr=nx_stds,fmt="o",label=f"{isetting}"+"x")
        # ax.errorbar(rotations_list,ny_means,yerr=ny_stds,fmt="<",label=f"{isetting}"+"y")
        if which_axis == "x":
            ax.errorbar(rotations_list,nx_means,yerr=nx_stds,fmt=">",label=f"g$_{{x}}$ {isetting}")
        elif which_axis == "y":
            ax.errorbar(rotations_list,ny_means,yerr=ny_stds,fmt=">",label=f"g$_{{y}}$ {isetting}")
        elif which_axis == "z":
            ax.errorbar(rotations_list,nz_means,yerr=nz_stds,fmt=">",label=f"g$_{{z}}$ {isetting}")
        else:
            ax.errorbar(rotations_list,nx_means,yerr=nx_stds,fmt=">",label=f"g$_{{x}}$ {isetting}")
            ax.errorbar(rotations_list,ny_means,yerr=ny_stds,fmt=">",label=f"g$_{{y}}$ {isetting}")
            # ax.errorbar(rotations_list,nz_means,yerr=nz_stds,fmt=">",label=f"g$_{{z}}$ {isetting}")

    # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
    # ax.set_ylabel(r" g_{x}", fontsize=22)
    # ax.set_ylabel(r" g_{y}", fontsize=22)
    ax.set_ylabel(f" g$_{{{which_axis}}}$ [ms$^{{-2}}$]", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    # ax.set_aspect('equal')
    ax.legend()
    # ax.set_xlim(0,360)
    # ax.set_ylim(0,10)
    # ax.set_yscale("log")
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/accelerometerxyz_finen{which_axis}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/accelerometerxyz_finen{which_axis}.pdf",transparent=False,bbox_inches='tight')
    plt.close()
# plot_accelerometer_variation_xyz(df,rotations_list,which_axis="x")
# plot_accelerometer_variation_xyz(df,rotations_list,which_axis="y")
# plot_accelerometer_variation_xyz(df,rotations_list,which_axis="z")
# plot_accelerometer_variation_xyz(df,rotations_list,which_axis="xyz")