#!/usr/bin/env python3

import glob
import json
import re

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates

import logging

# logging.basicConfig(level=logging.INFO)
logging.basicConfig(level=logging.ERROR)

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']


from pathlib import Path
home = str(Path.home())

from util import tilt_angle, to_spherical


geo_to_IC0 = 1948.07 #distance from surface of IC geometry to IC z=0

from tilt_monitor_cavern import (get_icm_id_list, get_icm_id_list_monidaq,
get_depth_from_icm_id, get_mean_tilt_icm_id,
plot_tilt_time_icm_ids,plot_tilt_time_icm_ids_height_subplots, plot_mean_tilt_depth)

from geometry_util import (get_device_type, get_string_port_from_icm_id,
                            get_depth_from_icm_id, get_string_position_from_icm_id)

from rotation_monitor_cavern import (get_time_rotation)




def get_icm_id_from_string_port(string,port,string_map_list):
    icm_id = None
    string_map = [ifile for ifile in string_map_list if f"string_{string}_" in ifile][0]
    with open(string_map, 'r') as f:
        data = json.load(f)
        for device in data:
            if device["port"] == port:
                icm_id = device["icm_id"]
                continue
    return icm_id



def plot_BxByBz(string,device,icm_list,moni_files,plotFolder,geometry_list,string_map_list):
    fig = plt.figure(figsize=(8,5*len(icm_list)))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig)
    for i,i_icm in enumerate(icm_list[:]):
        ax = fig.add_subplot(gs[i,0])
        time_list, Bx_list, By_list, Bz_list, B_list, B_theta_list, B_phi_list = get_time_rotation(string,i_icm,moni_files,geometry_list)
        print(f"icm id {i_icm} has {len(time_list)} rotation measurements")
        time_list = pd.to_datetime(time_list, format='ISO8601')
        upgrade_string,port = get_string_port_from_icm_id(i_icm,string_map_list)
        depth = get_depth_from_icm_id(i_icm,geometry_list)
        position,upgrade_string2 = get_string_position_from_icm_id(i_icm,geometry_list)
        ax.text(0.05, 0.95, f"S{upgrade_string}:{port} OM {int(position):d} {depth:.0f} m", transform=ax.transAxes, ha='left', fontsize=5)
        ax.plot(time_list,Bx_list,'-o',markersize=2,label=f"{'Bx'}",alpha=0.5)
        ax.plot(time_list,By_list,'-o',markersize=2,label=f"{'By'}",alpha=0.5)
        ax.plot(time_list,Bz_list,'-o',markersize=2,label=f"{'Bz'}",alpha=0.5)
        # print(f"time list {time_list}")
        ax.tick_params(axis='both',which='both', direction='in', labelsize=5)
        ax.grid(True,alpha=0.6)
        # ax.set_ylim(0,2.6)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
        ax.set_ylabel(r" B$_{i}$ [$\mu$T]", fontsize=5)
        ax.set_xlabel(f"time [UTC]", fontsize=5)
        ax.tick_params("x", rotation=45) 
        # if i < len(icm_list)-1:
        #     plt.setp(ax.get_xticklabels(), visible=False)
        #     ax.set_xlabel(None) 
        ax.legend(fontsize=5)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/Bx_By_Bz{string}_{device}_height_large.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_BxBy(string,device,icm_list,moni_files,plotFolder,geometry_list,string_map_list):
    fig = plt.figure(figsize=(8,5*len(icm_list)))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig)
    for i,i_icm in enumerate(icm_list[:]):
        ax = fig.add_subplot(gs[i,0])
        time_list, Bx_list, By_list, Bz_list, B_list, B_theta_list, B_phi_list = get_time_rotation(string,i_icm,moni_files,geometry_list)
        print(f"icm id {i_icm} has {len(time_list)} rotation measurements")
        time_list = pd.to_datetime(time_list, format='ISO8601')
        upgrade_string,port = get_string_port_from_icm_id(i_icm,string_map_list)
        depth = get_depth_from_icm_id(i_icm,geometry_list)
        position,upgrade_string2 = get_string_position_from_icm_id(i_icm,geometry_list)
        ax.text(0.05, 0.95, f"S{upgrade_string}:{port} OM {int(position):d} {depth:.0f} m", transform=ax.transAxes, ha='left', fontsize=16)
        ax.plot(Bx_list,By_list,'-o',markersize=2,label=f"{'Bx'}",alpha=0.5)
        # print(f"time list {time_list}")
        ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
        ax.grid(True,alpha=0.6)
        # ax.set_ylim(0,2.6)
        # ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
        ax.set_xlabel(r" B$_{x}$ [$\mu$T]", fontsize=16)
        ax.set_ylabel(r" B$_{y}$ [$\mu$T]", fontsize=16)
        ax.tick_params("x", rotation=45) 
        # if i < len(icm_list)-1:
        #     plt.setp(ax.get_xticklabels(), visible=False)
        #     ax.set_xlabel(None) 
        ax.legend(fontsize=16)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/Bx_By{string}_{device}_height_large.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def plot_BxByBz3D(string,device,icm_list,moni_files,plotFolder,geometry_list,string_map_list):
    fig = plt.figure(figsize=(8,5*len(icm_list)))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    ax = fig.add_subplot(111, projection='3d')
    for i,i_icm in enumerate(icm_list[:]):
        # ax = fig.add_subplot(gs[i,0])
        time_list, Bx_list, By_list, Bz_list, B_list, B_theta_list, B_phi_list = get_time_rotation(string,i_icm,moni_files,geometry_list)
        print(f"icm id {i_icm} has {len(time_list)} rotation measurements")
        time_list = pd.to_datetime(time_list, format='ISO8601')
        upgrade_string,port = get_string_port_from_icm_id(i_icm,string_map_list)
        depth = get_depth_from_icm_id(i_icm,geometry_list)
        position,upgrade_string2 = get_string_position_from_icm_id(i_icm,geometry_list)
        ax.text(0.05, 0.05,0.95, f"S{upgrade_string}:{port} OM {int(position):d} {depth:.0f} m", transform=ax.transAxes, ha='left', fontsize=16)
        ax.scatter(Bx_list,By_list,Bz_list,label=f"S{upgrade_string}:{port} OM {int(position):d} {depth:.0f} m",alpha=0.5)
        # print(f"time list {time_list}")
        ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
        ax.grid(True,alpha=0.6)
        # ax.set_ylim(0,2.6)
        # ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
        ax.set_xlabel(r" B$_{x}$ [$\mu$T]", fontsize=16)
        ax.set_ylabel(r" B$_{y}$ [$\mu$T]", fontsize=16)
        ax.set_zlabel(r" B$_{z}$ [$\mu$T]", fontsize=16)
        ax.tick_params("x", rotation=45) 
        ax.set_xlim(-70,-40)
        ax.set_ylim(10,40)
        ax.set_zlim(45,75)
        # if i < len(icm_list)-1:
        #     plt.setp(ax.get_xticklabels(), visible=False)
        #     ax.set_xlabel(None) 
        ax.legend(fontsize=16)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/Bx_By_Bz_3D{string}_{device}_{str(icm_list[0])}_height_large.pdf",transparent=False,bbox_inches='tight')
    # plt.show()
    plt.close()

def plot_B(string,device,icm_list,moni_files,plotFolder,geometry_list,string_map_list):
    fig = plt.figure(figsize=(8,5*len(icm_list)))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig)
    for i,i_icm in enumerate(icm_list[:]):
        ax = fig.add_subplot(gs[i,0])
        time_list, Bx_list, By_list, Bz_list, B_list, B_theta_list, B_phi_list = get_time_rotation(string,i_icm,moni_files,geometry_list)
        print(f"icm id {i_icm} has {len(time_list)} rotation measurements")
        time_list = pd.to_datetime(time_list, format='ISO8601')
        upgrade_string,port = get_string_port_from_icm_id(i_icm,string_map_list)
        depth = get_depth_from_icm_id(i_icm,geometry_list)
        position,upgrade_string2 = get_string_position_from_icm_id(i_icm,geometry_list)
        ax.text(0.05, 0.95, f"S{upgrade_string}:{port} OM {int(position):d} {depth:.0f} m", transform=ax.transAxes, ha='left', fontsize=16)
        ax.plot(time_list,B_list,'-o',markersize=2,label=f"{'B'}",alpha=0.5)
        # print(f"time list {time_list}")
        ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
        ax.grid(True,alpha=0.6)
        # ax.set_ylim(0,2.6)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
        ax.set_ylabel(r" B [$\mu$T]", fontsize=16)
        ax.set_xlabel(f"time [UTC]", fontsize=16)
        ax.tick_params("x", rotation=45)
        ax.set_ylim(0,90) 
        # if i < len(icm_list)-1:
        #     plt.setp(ax.get_xticklabels(), visible=False)
        #     ax.set_xlabel(None) 
        ax.legend(fontsize=16)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/B{string}_{device}_height_large.pdf",transparent=False,bbox_inches='tight')
    plt.close()






def main() -> None:
    from sensor_measurement import xDOMSensorMeasurement
    plotFolder = home+"/research_ua/icecube/upgrade/upgrade_comissioning/plots/"
    moni_data_folder = "/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/moni-data/"
    moni_files = sorted(glob.glob(moni_data_folder+"/*.json"))
    monidaq_device_list_file = "/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/monidaq_device_list.json"
    string_map_dir = "/Users/epaudel/Library/CloudStorage/OneDrive-TheUniversityofAlabama/research_notes/notes/string_maps/"
    string_map_list = sorted(glob.glob(string_map_dir+"/*.json"))
    geometry_folder = "/Users/epaudel/research_ua/icecube/software/upgrade_commissioning_scripts/geometry" #upgrade commissioning scripts repo
    geometry_list = sorted(glob.glob(geometry_folder+"/string_*geometry.json"))
    # tilt_plot("92","mDOM",moni_files,plotFolder)
    # tilt_plot("92","DEgg",moni_files,plotFolder)
    #get list of all devices on fieldhub 92
    # icm_ids_92 = get_icm_id_list("92",["degg","mdom"],geometry_list)
    icm_ids_92 = get_icm_id_list("92",["mdom"],geometry_list)
    # icm_id_s92_5154 = [get_icm_id_from_string_port("92", "mDOM", 5154, geometry_list)]
    # plot_rotation_time_icm_ids_height_subplots_large("92","mDOM",icm_ids_92,moni_files,plotFolder,geometry_list,string_map_list)
    # plot_Bphi_time_icm_ids_height_subplots_large("92","mDOM",icm_ids_92,moni_files,plotFolder,geometry_list,string_map_list)
    # plot_B_time_icm_ids_height_subplots_large("92","mDOM",icm_ids_92,moni_files,plotFolder,geometry_list,string_map_list)
    icm_id_92_5154 = get_icm_id_from_string_port("92", 5154, string_map_list)
    print(f"icm id for string 92 mDOM port 5154: {icm_id_92_5154}")
    plot_BxBy("92","mDOM",[icm_id_92_5154],moni_files,plotFolder,geometry_list,string_map_list)
    # plot_BxByBz3D("92","mDOM",[icm_id_92_5154],moni_files,plotFolder,geometry_list,string_map_list)
    plot_B("92","mDOM",[icm_id_92_5154],moni_files,plotFolder,geometry_list,string_map_list)
if __name__ == "__main__":
    main()