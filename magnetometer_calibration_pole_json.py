#!/usr/bin/env python

import os
import glob

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import json
import numpy as np

from pathlib import Path
home = str(Path.home())

from calibrate_magnetometer_2d import corrected_ellipse

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = home+"/research_ua/icecube/upgrade/upgrade_comissioning/plots/"
# print(degg_list)

colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']


#measurement taken on August 25, 2025 for magnetometer yNorth --> 360 rotation
import pandas as pd


mdom = "mDOM_pole"
#at pole
B_horizontal = 16.8 #muT
Bvert = 51.6 #muT
B_total = 54.2 #muT



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





def plot_corrected_heading_360(bx_mean_list,by_mean_list,angles_list,label=""):
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
    ax.set_xlim(0,360)
    ax.set_ylim(0,360)
    ax.set_aspect('equal')
    ax.legend(ncols=2,fontsize=16)
    ax.set_xlabel(r" mDOM rotation [$^{\circ}$]", fontsize=20)
    ax.set_ylabel(r" $\phi$ [$^{\circ}$]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_Bxy_heading_Pole_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_Bxy_heading_Pole_{label}.pdf",transparent=False,bbox_inches='tight')
    plt.close()




def plot_B_calibrated(mag_x, mag_y, mag_z,angles_list,label=""):
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
    # ax.set_aspect('equal')
    ax.legend(loc="lower left",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_ylim(0,90)
    ax.set_ylabel(r" $B$ [$\mu$T]", fontsize=20)
    ax.set_xlabel(r" $\phi$ [$^{\circ}$]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}B_calibrated_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}B_calibrated_{label}.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def plot_Bxy_calibrated(mag_x, mag_y, mag_z,angles_list,label=""):
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
    ax.set_ylim(0,90)
    # ax.set_aspect('equal')
    ax.legend(loc="upper right",ncols=1,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_ylabel(r" $B_{xy}$ [$\mu$T]", fontsize=20)
    ax.set_xlabel(r" $\phi$ [$^{\circ}$]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_{mdom}BxBy_calibrated_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_{mdom}BxBy_calibrated_{label}.pdf",transparent=False,bbox_inches='tight')
    plt.close()





def plot_Bz(mag_z,angles_list,label=""):
    mag_z = np.array(mag_z)*10**6 #convert to microTesla
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    ax.hlines(Bvert,0,360,ls="--",lw=2.5,label=f"B$_{{z}}$ ({Bvert:.1f} ${{\u03bc}}$T)",alpha=1.0)
    ax.plot(angles_list,mag_z, "-o", c="b", label=f"{'raw'}", alpha=1)
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




def main() -> None:
    pole_rotation_measurement_file = home+"/research_ua/icecube/upgrade/upgrade_comissioning/scripts/pole_rotation_measurement_run4.json"

    with open(pole_rotation_measurement_file, "r") as f:
        b_map = json.load(f)

    bx_collect = np.concatenate((b_map["bx_PMT8"],b_map["bx_PMT8p"],b_map["bx_PMT9"],b_map["bx_PMT9p"],b_map["bx_PMT10"],b_map["bx_PMT10p"],b_map["bx_PMT11"],b_map["bx_PMT11p"],b_map["bx_PMT4"],b_map["bx_PMT4p"],b_map["bx_PMT5"],b_map["bx_PMT5p"],b_map["bx_PMT6"],b_map["bx_PMT6p"],b_map["bx_PMT7"],b_map["bx_PMT7p"],b_map["bx_PMT8"]))
    by_collect = np.concatenate((b_map["by_PMT8"],b_map["by_PMT8p"],b_map["by_PMT9"],b_map["by_PMT9p"],b_map["by_PMT10"],b_map["by_PMT10p"],b_map["by_PMT11"],b_map["by_PMT11p"],b_map["by_PMT4"],b_map["by_PMT4p"],b_map["by_PMT5"],b_map["by_PMT5p"],b_map["by_PMT6"],b_map["by_PMT6p"],b_map["by_PMT7"],b_map["by_PMT7p"],b_map["by_PMT8"]))
    bz_collect = np.concatenate((b_map["bz_PMT8"],b_map["bz_PMT8p"],b_map["bz_PMT9"],b_map["bz_PMT9p"],b_map["bz_PMT10"],b_map["bz_PMT10p"],b_map["bz_PMT11"],b_map["bz_PMT11p"],b_map["bz_PMT4"],b_map["bz_PMT4p"],b_map["bz_PMT5"],b_map["bz_PMT5p"],b_map["bz_PMT6"],b_map["bz_PMT6p"],b_map["bz_PMT7"],b_map["bz_PMT7p"],b_map["bz_PMT8"]))

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
    for istep in step_list:
        bx_mean = np.mean(b_map[f"bx_PMT{istep}"])
        by_mean = np.mean(b_map[f"by_PMT{istep}"])
        bz_mean = np.mean(b_map[f"bz_PMT{istep}"])
        bx_mean_list.append(bx_mean)
        by_mean_list.append(by_mean)
        bz_mean_list.append(bz_mean)

    plot_xy_calibrated(bx_mean_list[:], by_mean_list[:],label="coarse")
    plot_corrected_heading_360(bx_mean_list, by_mean_list,angles_list,label="fine")
    plot_corrected_heading_360(bx_mean_list, by_mean_list,angles_list,label="coarse")
    plot_B_calibrated(bx_mean_list, by_mean_list, bz_mean_list,angles_list,label="coarse")
    plot_Bxy_calibrated(bx_mean_list, by_mean_list, bz_mean_list,angles_list,label="coarse")
    plot_Bz(bz_mean_list,angles_list,label="coarse")

if __name__ == "__main__":
    main()