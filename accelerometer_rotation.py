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


drive = "/Users/epaudel/Library/CloudStorage/OneDrive-TheUniversityofAlabama/"
accelerometer_readings = drive+"research_notes/references/upgrade_comissioning/mDOM_orientations/accelerometer_measurement.xlsx"
import pandas as pd


# df = pd.read_excel(accelerometer_readings,header=0,usecols="A:Q")
# df_tilts = pd.read_excel(accelerometer_readings,header=0,usecols="R:V")
# rotations = df.columns.values
# tilts = df_tilts.columns.values
# # print(df.columns.values)

# def tilt_angle(nx,ny,nz):
#     '''
#     calculate tilt angle
#     '''
#     return np.arctan2(np.sqrt(nx*nx + ny*ny),nz)

# def to_spherical(x,y,z):
#     '''
#     converts cartesian coordinates (x,y,z) to spherical (r,theta,phi)
#     '''
#     r = np.sqrt(x*x+y*y+z*z)
#     theta = np.arctan2(np.sqrt(x*x + y*y),z)
#     phi = np.arctan2(y,x)
#     return r,theta,phi
# g_means = []
# g_stds = []
# gxy_means = []
# gxy_stds = []
# zenith_means = []
# zenith_stds = []
# azimuth_means = []
# azimuth_stds = []
# for irot in rotations[:]:
#     accelerometer_readings = df[[irot]].values    
#     # print("accelerometer reading",accelerometer_readings)
#     angles = []
#     zenith = []
#     azimuth = []
#     g = []
#     gxy = []
#     for imeas in accelerometer_readings[:]:
#         print(imeas[0])
#         # nx,ny,nz = imeas[0]
#         print(imeas[0].split(","))
#         nx,ny,nz = [float(jx.strip("[").strip("]")) for jx in imeas[0].split(",")]
#         # print(tilt_angle(nx,ny,nz)*180/np.pi)
#         angles.append(tilt_angle(nx,ny,nz)*180/np.pi)
#         zenith.append(to_spherical(nx,ny,nz)[1]*180/np.pi)
#         azimuth.append(to_spherical(nx,ny,nz)[2]*180/np.pi)
#         g.append(to_spherical(nx,ny,nz)[0]*10**6)
#         gxy.append(np.sqrt(nx**2+ny**2))
#     # print(f"Rotation {irot} mean {np.mean(angles):.2f} std {np.std(angles):.2f}")s
#     # print(f"Rotation {irot} mean zen {np.mean(zenith):.2f} std {np.std(zenith):.2f}"+
#     #       f" mean azi {np.mean(azimuth):.2f} std {np.std(azimuth):.2f}"+
#     #       f" mean B {np.mean(B):.2f} std {np.std(B):.2f}")
#     azimuth = [i+360 if i<0 else i for i in azimuth]
#     azimuth_means.append(np.mean(azimuth))
#     azimuth_stds.append(np.std(azimuth))
#     zenith_means.append(np.mean(zenith))
#     zenith_stds.append(np.std(azimuth))
#     g_means.append(np.mean(g))
#     g_stds.append(np.std(g))
#     gxy_means.append(np.mean(gxy))
#     gxy_stds.append(np.std(gxy))

# # print(df.iloc[:,:1].values)
# def plot_accelerometer_variation(B,g_stds,zenith,zenith_stds,azimuth,azimuth_stds):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
#     ax.errorbar(rotations,azimuth,yerr=azimuth_stds,fmt="o",label="azimuth")
#     # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
#     ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
#     ax.set_ylabel(r" azimuth of $\vec{g}$ [$^{\circ}$]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     ax.set_xticks(np.linspace(0,360,9))
#     ax.set_yticks(np.linspace(0,360,9))
#     ax.set_aspect('equal')
#     ax.legend()
#     ax.set_xlim(0,360)
#     ax.set_ylim(0,360)
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/accelerometerxyz_fine.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/accelerometerxyz_fine.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# # plot_accelerometer_variation(g_means,g_stds,zenith_means,zenith_stds,azimuth_means,azimuth_stds)



# def plot_gxy(gxy_means, Bx_stds):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
#     ax.errorbar(rotations,gxy_means,yerr=gxy_stds,fmt="o",label="azimuth")
#     # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
#     ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
#     ax.set_ylabel(r" gxy[$\mu$T$]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     # ax.set_xticks(np.linspace(0,360,9))
#     # ax.set_yticks(np.linspace(0,360,9))
#     # ax.set_aspect('equal')
#     ax.legend()
#     # ax.set_xlim(0,360)
#     # ax.set_ylim(0,360)
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/accelerometergxy_fine.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/accelerometergxy_fine.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# # plot_gxy(gxy_means,gxy_stds)

# print(tilts)
# print(df_tilts)
# g_means = []
# g_stds = []
# gxy_means = []
# gxy_stds = []
# zenith_means = []
# zenith_stds = []
# azimuth_means = []
# azimuth_stds = []
# for irot in tilts[:]:
#     accelerometer_readings = df_tilts[[irot]].values 
#     print(accelerometer_readings)   
#     # print("accelerometer reading",accelerometer_readings)
#     angles = []
#     zenith = []
#     azimuth = []
#     g = []
#     gxy = []
#     for imeas in accelerometer_readings[:]:
#         print(imeas[0])
#         # nx,ny,nz = imeas[0]
#         print(imeas[0].split(","))
#         nx,ny,nz = [float(jx.strip("[").strip("]")) for jx in imeas[0].split(",")]
#         # print(tilt_angle(nx,ny,nz)*180/np.pi)
#         angles.append(tilt_angle(nx,ny,nz)*180/np.pi)
#         zenith.append(to_spherical(nx,ny,nz)[1]*180/np.pi)
#         azimuth.append(to_spherical(nx,ny,nz)[2]*180/np.pi)
#         g.append(to_spherical(nx,ny,nz)[0]*10**6)
#         gxy.append(np.sqrt(nx**2+ny**2))
#     # print(f"Rotation {irot} mean {np.mean(angles):.2f} std {np.std(angles):.2f}")s
#     # print(f"Rotation {irot} mean zen {np.mean(zenith):.2f} std {np.std(zenith):.2f}"+
#     #       f" mean azi {np.mean(azimuth):.2f} std {np.std(azimuth):.2f}"+
#     #       f" mean B {np.mean(B):.2f} std {np.std(B):.2f}")
#     azimuth = [i+360 if i<0 else i for i in azimuth]
#     azimuth_means.append(np.mean(azimuth))
#     azimuth_stds.append(np.std(azimuth))
#     zenith_means.append(np.mean(zenith))
#     zenith_stds.append(np.std(azimuth))
#     g_means.append(np.mean(g))
#     g_stds.append(np.std(g))
#     gxy_means.append(np.mean(gxy))
#     gxy_stds.append(np.std(gxy))

# # print(df.iloc[:,:1].values)
# def plot_accelerometer_variation_tilts(B,g_stds,zenith,zenith_stds,azimuth,azimuth_stds):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
#     print(tilts)
#     print(azimuth)
#     # ax.errorbar(tilts,zenith,yerr=zenith_stds,fmt="o",label="zenith")
#     ax.errorbar(tilts,azimuth,yerr=azimuth_stds,fmt="o",label="azimuth")
#     # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
#     # ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
#     ax.set_xlabel(r" tilt", fontsize=22)
#     ax.set_ylabel(r" azimuth of $\vec{g}$ [$^{\circ}$]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     # ax.set_xticks(np.linspace(0,360,9))
#     # ax.set_yticks(np.linspace(0,360,9))
#     # ax.set_aspect('equal')
#     ax.legend()
#     # ax.set_xlim(0,360)
#     # ax.set_ylim(0,360)
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/accelerometerxyz_fine_tilts.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/accelerometerxyz_fine_tilts.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# plot_accelerometer_variation_tilts(g_means,g_stds,zenith_means,zenith_stds,azimuth_means,azimuth_stds)


df = pd.read_excel(accelerometer_readings,sheet_name='accelerometer_04152025',header=0,usecols="A:GC")#Aup arrow tab up by 44+-1 mm
df_magnetometer = pd.read_excel(accelerometer_readings,sheet_name='magnetometer_04152025',header=0,usecols="A:GC")#Aup arrow tab up by 44+-1 mm
# df_tilts = pd.read_excel(accelerometer_readings,header=0,usecols="R:V")
rotations = df.columns.values
# tilts = df_tilts.columns.values
# print(df.columns.values)

def get_column_names(name_string):
    return [iname for iname in df.columns.values if name_string in iname]

df_B_up = df[get_column_names("B_up")]
# print(df_B_up)
# print("B_up",get_column_names("B_up"))
rotations_list = [int(izen.split("_")[0]) for izen in df_B_up.columns.values]
# print([int(izen.split("_")[0]) for izen in df_B_up.columns.values])

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

def get_means(df):
    g_means = []
    g_stds = []
    gxy_means = []
    gxy_stds = []
    zenith_means = []
    zenith_stds = []
    azimuth_means = []
    azimuth_stds = []
    nx_means = []
    ny_means = []
    nz_means = []
    nx_stds = []
    ny_stds = []
    nz_stds = []
    for irot in df.columns.values[:]:
        # print(f"{irot} rotation")
        accelerometer_readings = df[[irot]].values    
        # print("accelerometer reading",accelerometer_readings)
        angles = []
        zenith = []
        azimuth = []
        g = []
        gxy = []
        nx_list = []
        ny_list = []
        nz_list = []
        for imeas in accelerometer_readings[:]:
            # print(imeas[0])
            # nx,ny,nz = imeas[0]
            # print(imeas[0].split(","))
            nx,ny,nz = [float(jx.strip("[").strip("]")) for jx in imeas[0].split(",")]
            # print(tilt_angle(nx,ny,nz)*180/np.pi)
            angles.append(tilt_angle(nx,ny,nz)*180/np.pi)
            zenith.append(to_spherical(nx,ny,nz)[1]*180/np.pi)
            azimuth.append(to_spherical(nx,ny,nz)[2]*180/np.pi)
            g.append(to_spherical(nx,ny,nz)[0]*10**6)
            gxy.append(np.sqrt(nx**2+ny**2))
            nx_list.append(nx)
            ny_list.append(ny)
            nz_list.append(nz)
            # print(f"azimuth {azimuth}")
        # print(f"Rotation {irot} mean {np.mean(angles):.2f} std {np.std(angles):.2f}")s
        # print(f"Rotation {irot} mean zen {np.mean(zenith):.2f} std {np.std(zenith):.2f}"+
        #       f" mean azi {np.mean(azimuth):.2f} std {np.std(azimuth):.2f}"+
        #       f" mean B {np.mean(B):.2f} std {np.std(B):.2f}")
        azimuth = [i+360 if i<0 else i for i in azimuth]
        azimuth_means.append(np.mean(azimuth))
        azimuth_stds.append(np.std(azimuth))
        zenith_means.append(np.mean(zenith))
        zenith_stds.append(np.std(azimuth))
        g_means.append(np.mean(g))
        g_stds.append(np.std(g))
        gxy_means.append(np.mean(gxy))
        gxy_stds.append(np.std(gxy))
        nx_means.append(np.mean(nx_list))
        ny_means.append(np.mean(ny_list))
        nz_means.append(np.mean(nz_list))
        nx_stds.append(np.std(nx_list))
        ny_stds.append(np.std(ny_list))
        nz_stds.append(np.std(nz_list))

    # print(azimuth_means,g_means)
    return nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds
nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(df_B_up)
# print(df.iloc[:,:1].values)
tab_labels = {"plane":"plane","A_up":"tab A","B_up":"tab B","C_up":"tab C","D_up":"tab D"}
def plot_accelerometer_variation(df,rotation_list,sensor):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    # print(len(rotations_list),len(azimuth),azimuth)
    for isetting in ["plane","A_up","B_up","C_up","D_up"][:1]:
        this_df = df[get_column_names(isetting)]
        _,_,_,_,_,_,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        ax.errorbar(rotations_list,azimuth_means,yerr=azimuth_stds,fmt="o",label=f"{tab_labels[isetting]}")
    # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
    ax.set_ylabel(r" azimuth angle of $\vec{g}$ [$^{\circ}$]", fontsize=22)
    ax.grid(True,alpha=0.6)
    ax.set_xticks(np.linspace(0,360,9))
    ax.set_yticks(np.linspace(0,360,9))
    # ax.set_aspect('equal')
    ax.legend(ncols=4)
    ax.set_xlim(0,360)
    ax.set_ylim(0,360)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/{sensor}xyz_fine.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/{sensor}xyz_fine.pdf",transparent=False,bbox_inches='tight')
    plt.close()
plot_accelerometer_variation(df,rotations_list,"accelerometer")
plot_accelerometer_variation(df_magnetometer,rotations_list,"magnetometer")

def plot_accelerometer_variation_zenith(df,rotation_list,sensor):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
    # print(len(rotations_list),len(azimuth),azimuth)
    for isetting in ["plane","A_up","B_up","C_up","D_up"][:1]:
        this_df = df[get_column_names(isetting)]
        _,_,_,_,_,_,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        print(f"#######################zenith means {min(zenith_means)} {max(zenith_means)} {np.mean(zenith_means)}")
        ax.errorbar(rotations_list,zenith_means,yerr=zenith_stds,fmt="o",label=f"{tab_labels[isetting]}")
    # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
    ax.set_ylabel(r" zenith angle of $\vec{g}$ [$^{\circ}$]", fontsize=22)
    ax.grid(True,alpha=0.6)
    ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    # ax.set_aspect('equal')
    ax.legend()
    ax.set_xlim(0,360)
    # ax.set_ylim(0,1.25)
    # ax.set_ylim(0,12)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/{sensor}xyz_fineZenith.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/{sensor}xyz_fineZenith.pdf",transparent=False,bbox_inches='tight')
    plt.close()
plot_accelerometer_variation_zenith(df,rotations_list,"accelerometer")
plot_accelerometer_variation_zenith(df_magnetometer,rotations_list,"magnetometer")



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
plot_gxy(gxy_means,gxy_stds)



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


def plot_ginclinations_mdoms(df):
    fig = plt.figure(figsize=(8*3,5))
    # print(time_stamps)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0])
    shift = 0
    xtick_positions = []
    xtick_labels = []
    for isetting in ["plane","A_up","B_up","C_up","D_up"][:1]:
        this_df = df[get_column_names(isetting)]
        nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        for izenith in zenith_means:
            mean_x = np.sin(izenith)
            mean_y = np.cos(izenith)
            mean_xy = np.sqrt(mean_x**2+mean_y**2)
            mean_x/=mean_xy
            mean_y/=mean_xy
            print(mean_x,mean_y)
            # ax.plot([0,mean_bx],[0,mean_by],"-",alpha=1)
            ax.arrow(0+shift,0,mean_x,mean_y,width=0.01,alpha=1)
            ax.plot(0+shift,0,"o",ms=8,alpha=1)
            xtick_positions.append(0+shift)
            xtick_labels.append(izenith)
            shift += 10
    
    # ax.hist(deltaT,histtype="step",lw=2.5,alpha=1)
    # ax.set_xticks(xtick_positions,xtick_labels)
    ax.set_yticks([])
    # ax.set_title(f"channel {ikey}",fontsize=22,loc="left")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)                        
    # ax.set_ylabel(r"count", fontsize=22)        
    # ax.set_ylim(3450,3650)
    #test select######
    # ax.set_xlim(0,4100)        
    # ax.set_ylim(0,5000)        
    ax.grid(True,alpha=0.6)
    # print(f"mask count {maskCount} {sum(maskCount)}")
    ax.set_xlabel(r"directions", fontsize=22)
    # ax.legend()
    plt.savefig(plotFolder+f"/LabSensorgInclinationmDOMs_fine.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/LabSensorgInclinationmDOMs_fine.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_ginclinations_mdoms(df)


def plot_binclinations_mdoms(df):
    fig = plt.figure(figsize=(8*3,5))
    # print(time_stamps)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0])
    shift = 0
    xtick_positions = []
    xtick_labels = []
    for isetting in ["plane","A_up","B_up","C_up","D_up"][:1]:
        this_df = df[get_column_names(isetting)]
        nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        for izenith in zenith_means:
            mean_x = np.sin(izenith)
            mean_y = np.cos(izenith)
            mean_xy = np.sqrt(mean_x**2+mean_y**2)
            mean_x/=mean_xy
            mean_y/=mean_xy
            print(mean_x,mean_y)
            # ax.plot([0,mean_bx],[0,mean_by],"-",alpha=1)
            ax.arrow(0+shift,0,mean_x,mean_y,width=0.01,alpha=1)
            ax.plot(0+shift,0,"o",ms=8,alpha=1)
            xtick_positions.append(0+shift)
            xtick_labels.append(izenith)
            shift += 10
    
    # ax.hist(deltaT,histtype="step",lw=2.5,alpha=1)
    ax.set_xticks(xtick_positions,xtick_labels)
    ax.set_yticks([])
    # ax.set_title(f"channel {ikey}",fontsize=22,loc="left")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)                        
    # ax.set_ylabel(r"count", fontsize=22)        
    # ax.set_ylim(3450,3650)
    #test select######
    # ax.set_xlim(0,4100)        
    # ax.set_ylim(0,5000)        
    ax.grid(True,alpha=0.6)
    # print(f"mask count {maskCount} {sum(maskCount)}")
    ax.set_xlabel(r"directions", fontsize=22)
    # ax.legend()
    plt.savefig(plotFolder+f"/LabSensorbInclinationmDOMs_fine.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/LabSensorbInclinationmDOMs_fine.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_binclinations_mdoms(df_magnetometer)

def plot_gorientations_mdoms(df):
    fig = plt.figure(figsize=(8,5))
    # print(time_stamps)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0])
    shift = 0
    xtick_positions = []
    xtick_labels = []
    for isetting in ["plane","A_up","B_up","C_up","D_up"][:1]:
        this_df = df[get_column_names(isetting)]
        nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        for igx,igy in zip(nx_means,ny_means):
            mean_bx = np.mean(igx)
            mean_by = np.mean(igy)
            mean_bxby = np.sqrt(mean_bx**2+mean_by**2)
            mean_bx/=mean_bxby
            mean_by/=mean_bxby
            # ax.plot([0,mean_bx],[0,mean_by],"-",alpha=1)
            ax.arrow(0+shift,0,mean_bx,mean_by,width=0.01,alpha=1)
            ax.plot(0+shift,0,"o",ms=8,alpha=1)
            xtick_positions.append(0+shift)
            # xtick_labels.append(isetting)
            shift += 10
    
    # ax.hist(deltaT,histtype="step",lw=2.5,alpha=1)
    # ax.set_xticks(xtick_positions,xtick_labels)
    ax.set_yticks([])
    # ax.set_title(f"channel {ikey}",fontsize=22,loc="left")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)                        
    # ax.set_ylabel(r"count", fontsize=22)        
    # ax.set_ylim(3450,3650)
    #test select######
    # ax.set_xlim(0,4100)        
    # ax.set_ylim(0,5000)        
    ax.grid(True,alpha=0.6)
    # print(f"mask count {maskCount} {sum(maskCount)}")
    ax.set_xlabel(r"directions", fontsize=22)
    # ax.legend()
    plt.savefig(plotFolder+f"/LabSensorgOrientationsmDOMs_fine.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/LabSensorgOrientationsmDOMs_fine.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_gorientations_mdoms(df)


def plot_borientations_fine_mdoms(df):
    fig = plt.figure(figsize=(8,5))
    # print(time_stamps)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0])
    shift = 0
    xtick_positions = []
    xtick_labels = []
    for isetting in ["plane","A_up","B_up","C_up","D_up"][:1]:
        this_df = df[get_column_names(isetting)]
        nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
        # print(np.linspace(0,360,37),np.zeros(len(nx_means)))
        # print(nx_means,ny_means)
        # ax.quiver([np.linspace(0,360,37),np.zeros(len(nx_means))],[nx_means,ny_means],angles="uv")
        
        print(f"azimuth{azimuth_means[:]}")
        for igx,igy in zip(nx_means[:],ny_means[:]):
            print(f"igx {igx}")
            mean_bx = igx
            mean_by = igy
            mean_bxby = np.sqrt(mean_bx**2+mean_by**2)
            mean_bx/=mean_bxby*10
            mean_by/=mean_bxby*10
            print(f"means {mean_bx} {mean_by}")
            # ax.plot([0,mean_bx],[0,mean_by],"-",alpha=1)
            # ax.arrow(0+shift,0,mean_bx+shift,mean_by,width=0.01,alpha=1)
            ax.annotate("", xytext=(0+shift, 0), xy=(mean_bx+shift, mean_by),arrowprops=dict(arrowstyle="->"))
            ax.plot(0+shift,0,"o",ms=8,alpha=1)
            xtick_positions.append(0+shift)
            # xtick_labels.append(isensor)
            shift += 10
    
    
    # ax.hist(deltaT,histtype="step",lw=2.5,alpha=1)
    # ax.set_xticks(xtick_positions,xtick_labels)
    ax.set_yticks([])
    # ax.set_title(f"channel {ikey}",fontsize=22,loc="left")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)                        
    # ax.set_ylabel(r"count", fontsize=22)        
    # ax.set_ylim(3450,3650)
    #test select######
    # ax.set_xlim(0,360)        
    # ax.set_ylim(0,5000)        
    ax.grid(True,alpha=0.6)
    # print(f"mask count {maskCount} {sum(maskCount)}")
    ax.set_xlabel(r"directions", fontsize=22)
    # ax.legend()
    plt.savefig(plotFolder+f"/LabSensorbOrientationsmDOMs_fine.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/LabSensorbOrientationsmDOMs_fine.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_borientations_fine_mdoms(df_magnetometer)













################mini mainboard##############
############################################
############################################
# df = pd.read_excel(accelerometer_readings,sheet_name='accelerometer_05132025MiniMB_re',header=0,usecols="A:AK")#Aup arrow tab up by 44+-1 mm
# # df = pd.read_excel(accelerometer_readings,sheet_name='accelerometer_05132025MiniMBDup',header=0,usecols="A:AK")#Aup arrow tab up by 44+-1 mm
# # df_tilts = pd.read_excel(accelerometer_readings,header=0,usecols="R:V")
# rotations = df.columns.values
# # tilts = df_tilts.columns.values
# # print(df.columns.values)

# def get_column_names(name_string):
#     return [iname for iname in df.columns.values if name_string in iname]

# # df_B_up = df[get_column_names("")]
# # print(df_B_up)
# # print("B_up",get_column_names("B_up"))
# # rotations_list = [int(izen.split("_")[0]) for izen in df_B_up.columns.values]
# # print([int(izen.split("_")[0]) for izen in df_B_up.columns.values])

# def tilt_angle(nx,ny,nz):
#     '''
#     calculate tilt angle
#     '''
#     return np.arctan2(np.sqrt(nx*nx + ny*ny),nz)

# def to_spherical(x,y,z):
#     '''
#     converts cartesian coordinates (x,y,z) to spherical (r,theta,phi)
#     '''
#     r = np.sqrt(x*x+y*y+z*z)
#     theta = np.arctan2(np.sqrt(x*x + y*y),z)
#     phi = np.arctan2(y,x)
#     return r,theta,phi

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



# def plot_gxy(gxy_means, Bx_stds):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
#     # print(len(rotations_list),len(azimuth),azimuth)
#     for isetting in ["plane","A_up","B_up","C_up","D_up"]:
#         this_df = df[get_column_names(isetting)]
#         nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
#         ax.errorbar(rotations_list,gxy_means,yerr=gxy_stds,fmt="o",label=f"{isetting}")
#     # ax.errorbar(rotations,gxy_means,yerr=gxy_stds,fmt="o",label="azimuth")
#     # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
#     ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
#     ax.set_ylabel(r" gxy[$\mu$T$]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     # ax.set_xticks(np.linspace(0,360,9))
#     # ax.set_yticks(np.linspace(0,360,9))
#     # ax.set_aspect('equal')
#     ax.legend()
#     # ax.set_xlim(0,360)
#     # ax.set_ylim(0,360)
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/accelerometergxy_fine.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/accelerometergxy_fine.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# # plot_gxy(gxy_means,gxy_stds)



# def plot_accelerometer_variation_xyz(df,rotation_list):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
#     # print(len(rotations_list),len(azimuth),azimuth)
#     for isetting in ["plane","A_up","B_up","C_up","D_up"][:1]:
#         this_df = df[get_column_names(isetting)]
#         nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
#         # ax.errorbar(rotations_list,nx_means,yerr=nx_stds,fmt="o",label=f"{isetting}"+"x")
#         # ax.errorbar(rotations_list,ny_means,yerr=ny_stds,fmt="<",label=f"{isetting}"+"y")
#         ax.errorbar(rotations_list,nz_means,yerr=nz_stds,fmt=">",label=f"g$_{{z}}$ {isetting}")
#     # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
#     ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
#     # ax.set_ylabel(r" g_{x}", fontsize=22)
#     # ax.set_ylabel(r" g_{y}", fontsize=22)
#     ax.set_ylabel(r" g$_{z}$ [ms$^{-2}$]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     # ax.set_xticks(np.linspace(0,360,9))
#     # ax.set_yticks(np.linspace(0,360,9))
#     # ax.set_aspect('equal')
#     ax.legend()
#     # ax.set_xlim(0,360)
#     # ax.set_ylim(0,10)
#     # ax.set_yscale("log")
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/accelerometerxyz_finenxnynz.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/accelerometerxyz_finenxnynz.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# # plot_accelerometer_variation_xyz(df,rotations_list)
# axis_label = {"x"}
# def plot_accelerometer_variation_xyz(df,rotation_list,which_axis):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     # ax.errorbar(rotations,zenith,yerr=zenith_stds,fmt="o",label="zenith")
#     # print(len(rotations_list),len(azimuth),azimuth)
#     for isetting in ["plane","A_up","B_up","C_up","D_up"][1:2]:
#         this_df = df[get_column_names(isetting)]
#         nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds = get_means(this_df)
#         # ax.errorbar(rotations_list,nx_means,yerr=nx_stds,fmt="o",label=f"{isetting}"+"x")
#         # ax.errorbar(rotations_list,ny_means,yerr=ny_stds,fmt="<",label=f"{isetting}"+"y")
#         if which_axis == "x":
#             ax.errorbar(rotations_list,nx_means,yerr=nx_stds,fmt=">",label=f"g$_{{x}}$ {isetting}")
#         elif which_axis == "y":
#             ax.errorbar(rotations_list,ny_means,yerr=ny_stds,fmt=">",label=f"g$_{{y}}$ {isetting}")
#         elif which_axis == "z":
#             ax.errorbar(rotations_list,nz_means,yerr=nz_stds,fmt=">",label=f"g$_{{z}}$ {isetting}")
#         else:
#             ax.errorbar(rotations_list,nx_means,yerr=nx_stds,fmt=">",label=f"g$_{{x}}$ {isetting}")
#             ax.errorbar(rotations_list,ny_means,yerr=ny_stds,fmt=">",label=f"g$_{{y}}$ {isetting}")
#             # ax.errorbar(rotations_list,nz_means,yerr=nz_stds,fmt=">",label=f"g$_{{z}}$ {isetting}")

#     # ax.errorbar(rotations,B,yerr=g_stds,fmt="o",label="B")
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
#     ax.set_xlabel(r" rotation [$^{\circ}$]", fontsize=22)
#     # ax.set_ylabel(r" g_{x}", fontsize=22)
#     # ax.set_ylabel(r" g_{y}", fontsize=22)
#     ax.set_ylabel(f" g$_{{{which_axis}}}$ [ms$^{{-2}}$]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     # ax.set_xticks(np.linspace(0,360,9))
#     # ax.set_yticks(np.linspace(0,360,9))
#     # ax.set_aspect('equal')
#     ax.legend()
#     # ax.set_xlim(0,360)
#     # ax.set_ylim(0,10)
#     # ax.set_yscale("log")
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/accelerometerxyz_finen{which_axis}.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/accelerometerxyz_finen{which_axis}.pdf",transparent=False,bbox_inches='tight')
#     plt.close()
# # plot_accelerometer_variation_xyz(df,rotations_list,which_axis="x")
# # plot_accelerometer_variation_xyz(df,rotations_list,which_axis="y")
# # plot_accelerometer_variation_xyz(df,rotations_list,which_axis="z")
# # plot_accelerometer_variation_xyz(df,rotations_list,which_axis="xyz")