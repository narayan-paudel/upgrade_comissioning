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
from geometry_util import (get_device_type, get_string_port_from_icm_id,
                            get_depth_from_icm_id, get_string_position_from_icm_id,
                            get_icm_id_from_string_port)


geo_to_IC0 = 1948.07 #distance from surface of IC geometry to IC z=0







def plot_orientation_moni(orientation_list,time_list,board_name,hostname,port,icm_id):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    bins = np.linspace(8.8,10.1,40)
    time_list = pd.to_datetime(time_list, format='ISO8601')
    ax.plot(time_list,orientation_list,'o',markersize=2,alpha=0.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    ax.grid(True,alpha=0.6)
    # ax.legend(,ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    # ax.set_ylim(-180,180)
    # ax.set_ylim(0,360)
    # ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    # ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.set_ylabel(f" orientation angle [\u00b0]", fontsize=20)
    ax.set_xlabel(f"time [UTC]", fontsize=20)
    ax.tick_params("x", rotation=45) 
    # plt.savefig(plotFolder+f"/orientation_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_vs_time_{board_name}_{hostname}_{port}_{icm_id}.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def orientation_time_from_moni(moni_files):
    for ifile in moni_files[:]:
        with open(ifile, "r") as fh:
            moni_data = json.load(fh)
            # print(f"moni data {moni_data}")
            count = 0
            tilt_list = []
            orientation_list = []
            time_list = []
            for ievent in moni_data[:]:
                # print(f"{ievent['moni_common']['accelerometer']}")
                gx,gy,gz = ievent['moni_common']['accelerometer']
                bx,by,bz = ievent['moni_common']['magnetometer']
                g,g_theta,g_phi = to_spherical(gx,gy,gz)
                B,B_theta,B_phi = to_spherical(bx,by,bz)
                g_tilt = tilt_angle(gx,gy,gz)
                # print(f"tilt {np.rad2deg(g_tilt)} degrees")
                # print(f"azimuth {np.rad2deg(g_phi)} degrees")
                time = ievent["header"]['utc_time']
                board_name = ievent['header']['board_name']
                hostname = ievent['header']['hostname']
                port = ievent['header']['port']
                icm_id = ievent['header']['icm_id']
                # print(f" board name{ievent['header']['board_name']} hostname {ievent['header']['hostname']} port {ievent['header']['port']} icm id {ievent['header']['icm_id']} ")
                tilt_list.append(np.rad2deg(g_tilt))
                orientation_list.append(np.rad2deg(B_phi))
                time_list.append(time)
    return tilt_list, orientation_list, time_list

def plot_tilt_moni(tilt_list,time_list,board_name,hostname,port,icm_id,plotFolder):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    bins = np.linspace(8.8,10.1,40)
    time_list = pd.to_datetime(time_list, format='ISO8601')
    ax.plot(time_list,tilt_list,'o',markersize=2,alpha=0.5)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    ax.grid(True,alpha=0.6)
    # ax.legend(,ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    # ax.set_ylim(0.5,0.7)
    # ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    # print(ax.get_xticklabels(minor=False))
    # labels = [label.get_text()[0:10] for label in ax.get_xticklabels(minor=False)]
    # location = [label.get_position() for label in ax.get_xticklabels(minor=False)]
    # label_loc = [text.Text(iloc,0,ilabel) for iloc, ilabel in enumerate(zip(location,labels))]
    # ax.set_xticklabels(label_loc,rotation=45, ha="right")
    ###################
    # ax.set_xticklabels(labels,rotation=45, ha="right")
    ################
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    # ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    # print(f"xtick labels {labels}, {location}")
    ax.set_ylabel(f" tilt angle [\u00b0]", fontsize=20)
    ax.set_xlabel(f"time [UTC]", fontsize=20)
    ax.tick_params("x", rotation=45) 
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_orientation_moni(orientation_list,time_list,board_name,hostname,port,icm_id,plotFolder):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    bins = np.linspace(8.8,10.1,40)
    time_list = pd.to_datetime(time_list, format='ISO8601')
    ax.plot(time_list,orientation_list,'o',markersize=2,alpha=0.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    ax.grid(True,alpha=0.6)
    # ax.legend(,ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    # ax.set_ylim(-180,180)
    # ax.set_ylim(0,360)
    # ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    # ax.xaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.set_ylabel(f" orientation angle [\u00b0]", fontsize=20)
    ax.set_xlabel(f"time [UTC]", fontsize=20)
    ax.tick_params("x", rotation=45) 
    plt.savefig(plotFolder+f"/orientation_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_vs_time_{board_name}_{hostname}_{port}_{icm_id}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def tilt_plot(string,device,moni_files,plotFolder):
    tilt_list = []
    time_list = []
    for ifile in moni_files[:]:
        if f"fieldhub{string}_{device}" in ifile:
            with open(ifile, "r") as fh:
                moni_data = json.load(fh)
                # print(f"moni data {moni_data}")
                count = 0
                for ievent in moni_data[:]:
                    time = ievent["header"]['utc_time']
                    board_name = ievent['header']['board_name']
                    hostname = ievent['header']['hostname']
                    port = ievent['header']['port']
                    icm_id = ievent['header']['icm_id']
                    # print(f"{ievent['moni_common']['accelerometer']}")
                    gx,gy,gz = ievent['moni_common']['accelerometer']
                    bx,by,bz = ievent['moni_common']['magnetometer']
                    g,g_theta,g_phi = to_spherical(gx,gy,gz)
                    B,B_theta,B_phi = to_spherical(bx,by,bz)
                    g_tilt = tilt_angle(gx,gy,gz)
                    # if g_tilt > np.deg2rad(2.1):
                    #     print(f"Large tilt found in file {board_name} {hostname} {port} {icm_id}: {np.rad2deg(g_tilt)} degrees")
                    # print(f"tilt {np.rad2deg(g_tilt)} degrees")
                    # print(f"azimuth {np.rad2deg(g_phi)} degrees")
                    # print(f" board name{ievent['header']['board_name']} hostname {ievent['header']['hostname']} port {ievent['header']['port']} icm id {ievent['header']['icm_id']} ")
                    if int(port) ==5154:
                        tilt_list.append(np.rad2deg(g_tilt))
                        # orientation_list.append(np.rad2deg(B_phi))
                        time_list.append(time)
                # print(f"plotting file {string} {device}")
    plot_tilt_moni(tilt_list,time_list,device,f"fieldhub{string}",f"all",f"all",plotFolder)
    # plot_orientation_moni(orientation_list,time_list,board_name,hostname,port,icm_id,plotFolder)

# tilt_plot("87","mDOM")
# tilt_plot("88","mDOM")
# tilt_plot("89","mDOM")
# tilt_plot("90","mDOM")
# tilt_plot("91","mDOM")
# tilt_plot("92","mDOM")




def orientation_plot(string,device):
    orientation_list = []
    time_list = []
    for ifile in moni_files[:]:
        if f"fieldhub{string}_{device}" in ifile:
            with open(ifile, "r") as fh:
                moni_data = json.load(fh)
                # print(f"moni data {moni_data}")
                count = 0
                for ievent in moni_data[:]:
                    time = ievent["header"]['utc_time']
                    board_name = ievent['header']['board_name']
                    hostname = ievent['header']['hostname']
                    port = ievent['header']['port']
                    icm_id = ievent['header']['icm_id']
                    # print(f"{ievent['moni_common']['accelerometer']}")
                    gx,gy,gz = ievent['moni_common']['accelerometer']
                    bx,by,bz = ievent['moni_common']['magnetometer']
                    g,g_theta,g_phi = to_spherical(gx,gy,gz)
                    B,B_theta,B_phi = to_spherical(bx,by,bz)
                    g_tilt = tilt_angle(gx,gy,gz)
                    # print(f"tilt {np.rad2deg(g_tilt)} degrees")
                    # print(f"azimuth {np.rad2deg(g_phi)} degrees")
                    # print(f" board name{ievent['header']['board_name']} hostname {ievent['header']['hostname']} port {ievent['header']['port']} icm id {ievent['header']['icm_id']} ")
                    orientation_list.append(np.rad2deg(B_phi))
                    time_list.append(time)
                # print(f"plotting file {string} {device}")
    plot_orientation_moni(orientation_list,time_list,device,f"fieldhub{string}",f"all",f"all",plotFolder)


def get_icm_id_list(string,device_list,geometry_list):
    icm_id_list = []
    string_geometry_file = [igeometry for igeometry in geometry_list if f"string_{string}_" in igeometry][0]
    with open(string_geometry_file, 'r') as f:
        geometry_data = json.load(f)
        for device in geometry_data[0]['devices']:
            if device['type'] in device_list:
                icm_id_list.append(device['icm_id'])
    return icm_id_list


def get_icm_id_list_depth(string,device_list,geometry_list):
    icm_id_list = []
    string_geometry_file = [igeometry for igeometry in geometry_list if f"string_{string}_" in igeometry][0]
    with open(string_geometry_file, 'r') as f:
        geometry_data = json.load(f)
        for device in geometry_data[0]['devices']:
            if device['type'] in device_list:
                icm_id_list.append(device['icm_id'])
    return icm_id_list

def get_icm_id_list_monidaq(string,select_devices,monidaq_device_list):
    icm_id_list = []
    with open(monidaq_device_list, "r") as fh:
        device_list = json.load(fh)
    for ifieldhub in device_list[:]:
        # print(f"fieldhub id {ifieldhub['fieldhub_id']}")
        if str(string) in ifieldhub["fieldhub_id"]:
            # print(f"found fieldhub {ifieldhub['fieldhub_id']} with {len(ifieldhub['devices'])} devices")
            for idevice in ifieldhub["devices"]:
                # print(f"device {idevice['board_name'].lower()} with icm id {idevice['icm_id']}")
                if idevice["board_name"].lower() in select_devices:
                    icm_id_list.append(idevice["icm_id"])
    return icm_id_list


def tilt_plot_port(string,device,port_list,moni_files,plotFolder):
    tilt_list = []
    time_list = []
    for ifile in moni_files[:]:
        if f"fieldhub{string}_{device}" in ifile:
            with open(ifile, "r") as fh:
                moni_data = json.load(fh)
                # print(f"moni data {moni_data}")
                count = 0
                for ievent in moni_data[:]:
                    time = ievent["header"]['utc_time']
                    board_name = ievent['header']['board_name']
                    hostname = ievent['header']['hostname']
                    port = ievent['header']['port']
                    icm_id = ievent['header']['icm_id']
                    if int(port) not in port_list:
                        continue
                    # print(f"{ievent['moni_common']['accelerometer']}")
                    gx,gy,gz = ievent['moni_common']['accelerometer']
                    bx,by,bz = ievent['moni_common']['magnetometer']
                    g,g_theta,g_phi = to_spherical(gx,gy,gz)
                    B,B_theta,B_phi = to_spherical(bx,by,bz)
                    g_tilt = tilt_angle(gx,gy,gz)
                    # if g_tilt > np.deg2rad(2.1):
                    #     print(f"Large tilt found in file {board_name} {hostname} {port} {icm_id}: {np.rad2deg(g_tilt)} degrees")
                    # print(f"tilt {np.rad2deg(g_tilt)} degrees")
                    # print(f"azimuth {np.rad2deg(g_phi)} degrees")
                    # print(f" board name{ievent['header']['board_name']} hostname {ievent['header']['hostname']} port {ievent['header']['port']} icm id {ievent['header']['icm_id']} ")
                    tilt_list.append(np.rad2deg(g_tilt))
                    # orientation_list.append(np.rad2deg(B_phi))
                    time_list.append(time)
                # print(f"plotting file {string} {device}")
    plot_tilt_moni(tilt_list,time_list,device,f"fieldhub{string}",string,str(port_list[0]),plotFolder)
    # plot_orientation_moni(orientation_list,time_list,board_name,hostname,port,icm_id,plotFolder)



def tilt_plot_port_range(string,device,low_port,high_port,moni_files,plotFolder):
    tilt_list = []
    time_list = []
    for ifile in moni_files[:]:
        if f"fieldhub{string}_{device}" in ifile:
            with open(ifile, "r") as fh:
                moni_data = json.load(fh)
                # print(f"moni data {moni_data}")
                count = 0
                for ievent in moni_data[:]:
                    time = ievent["header"]['utc_time']
                    board_name = ievent['header']['board_name']
                    hostname = ievent['header']['hostname']
                    port = ievent['header']['port']
                    icm_id = ievent['header']['icm_id']
                    if int(port) < low_port or int(port) > high_port:
                        continue
                    # print(f"{ievent['moni_common']['accelerometer']}")
                    gx,gy,gz = ievent['moni_common']['accelerometer']
                    bx,by,bz = ievent['moni_common']['magnetometer']
                    g,g_theta,g_phi = to_spherical(gx,gy,gz)
                    B,B_theta,B_phi = to_spherical(bx,by,bz)
                    g_tilt = tilt_angle(gx,gy,gz)
                    # if g_tilt > np.deg2rad(2.1):
                    #     print(f"Large tilt found in file {board_name} {hostname} {port} {icm_id}: {np.rad2deg(g_tilt)} degrees")
                    # print(f"tilt {np.rad2deg(g_tilt)} degrees")
                    # print(f"azimuth {np.rad2deg(g_phi)} degrees")
                    # print(f" board name{ievent['header']['board_name']} hostname {ievent['header']['hostname']} port {ievent['header']['port']} icm id {ievent['header']['icm_id']} ")
                    tilt_list.append(np.rad2deg(g_tilt))
                    # orientation_list.append(np.rad2deg(B_phi))
                    time_list.append(time)
                # print(f"plotting file {string} {device}")
    plot_tilt_moni(tilt_list,time_list,device,f"fieldhub{string}",string,str(low_port)+"_"+str(high_port),plotFolder)








def get_time_tilt(string,icm_id,moni_files,geometry_list):
    tilt_list = []
    time_list = []
    device_type = get_device_type(icm_id,geometry_list)
    logging.info(f"getting tilt for icm id {icm_id} of type {device_type} on string {string}")
    
    for ifile in moni_files[:]:
        if f"fieldhub{string}_{device_type}_{icm_id}" in ifile:
            with open(ifile, "r") as fh:
                moni_data = json.load(fh)
                # print(f"moni data {moni_data}")
                count = 0
                for ievent in moni_data[:]:
                    time = ievent["header"]['utc_time']
                    board_name = ievent['header']['board_name']
                    hostname = ievent['header']['hostname']
                    port = ievent['header']['port']
                    icm_id_file = ievent['header']['icm_id']
                    if icm_id_file != icm_id:
                        continue
                    # print(f"{ievent['moni_common']['accelerometer']}")
                    gx,gy,gz = ievent['moni_common']['accelerometer']
                    bx,by,bz = ievent['moni_common']['magnetometer']
                    g,g_theta,g_phi = to_spherical(gx,gy,gz)
                    B,B_theta,B_phi = to_spherical(bx,by,bz)
                    g_tilt = tilt_angle(gx,gy,gz)
                    # if g_tilt > np.deg2rad(2.1):
                    #     print(f"Large tilt found in file {board_name} {hostname} {port} {icm_id}: {np.rad2deg(g_tilt)} degrees")
                    # print(f"tilt {np.rad2deg(g_tilt)} degrees")
                    # print(f"azimuth {np.rad2deg(g_phi)} degrees")
                    # print(f" board name{ievent['header']['board_name']} hostname {ievent['header']['hostname']} port {ievent['header']['port']} icm id {ievent['header']['icm_id']} ")
                    tilt_list.append(np.rad2deg(g_tilt))
                    # orientation_list.append(np.rad2deg(B_phi))
                    time_list.append(time)
                # print(f"plotting file {string} {device}")
    return time_list, tilt_list

def get_mean_tilt_icm_id(string,icm_id,moni_files,geometry_list,date_start=None,date_end=None):
    time_list, tilt_list = get_time_tilt(string,icm_id,moni_files,geometry_list)
    time_list = pd.to_datetime(time_list, format='ISO8601')
    tilt_array = np.array(tilt_list)
    print(f"ICM ID {icm_id} has {tilt_array} measurements")
    if date_start is not None:
        date_start = pd.to_datetime(date_start,format="ISO8601")
        mask_start = time_list >= date_start
        # print(f"mask start {mask_start}")
        time_list = time_list[mask_start]
        tilt_array = tilt_array[mask_start]
    if date_end is not None:
        date_end = pd.to_datetime(date_end,format="ISO8601")
        mask_end = time_list <= date_end
        # print(f"mask end {mask_end}")
        time_list = time_list[mask_end]
        tilt_array = tilt_array[mask_end]
    if len(tilt_array) > 0:
        print(f" min max time {min(time_list)} {max(time_list)}")
        # print(f"mean tilt for icm id {icm_id} on string {string} between {date_start} and {date_end} has {tilt_array} measurements")
        mean_tilt = np.mean(tilt_array)
        std_tilt = np.std(tilt_array)
        return mean_tilt,std_tilt
    else:
        print(f"No measurements found for ICM ID {icm_id} on string {string} between {date_start} and {date_end}")
        return np.nan,np.nan

def plot_mean_tilt_depth(string,icm_list,moni_files,plotFolder,geometry_list,date_start=None,date_end=None):
    depth_list = []
    mean_tilt_list = []
    std_tilt_list = []
    for i_icm in icm_list[:]:
        depth = get_depth_from_icm_id(i_icm,geometry_list)
        mean_tilt,std_tilt = get_mean_tilt_icm_id(string,i_icm,moni_files,geometry_list,date_start=date_start,date_end=date_end)
        print(f"icm id {i_icm} has mean tilt {mean_tilt} degrees at depth {depth} m")
        depth_list.append(depth)
        mean_tilt_list.append(mean_tilt)
        std_tilt_list.append(std_tilt)
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    ax.plot(depth_list,mean_tilt_list,'o',markersize=6,label=f"{"tilt"}",alpha=0.5)
    # ax.errorbar(depth_list,mean_tilt_list,yerr=std_tilt_list,fmt='o',markersize=6,label=f"{"tilt"}",alpha=0.5)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.set_ylabel(f" mean tilt angle [\u00b0]", fontsize=20)
    ax.set_xlabel(f" depth [m]", fontsize=20)
    ax.legend(fontsize=16)
    plt.savefig(plotFolder+f"/mean_tilt_vs_depth_{string}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_tilt_time_icm_ids(string,device,icm_list,moni_files,plotFolder,geometry_list):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    for i_icm in icm_list[:]:
        time_list, tilt_list = get_time_tilt(string,device,i_icm,moni_files,geometry_list)
        print(f"icm id {i_icm} has {len(time_list)} tilt measurements")
        time_list = pd.to_datetime(time_list, format='ISO8601')
        ax.plot(time_list,tilt_list,'-o',markersize=2,label=f"ICM {i_icm}",alpha=0.5)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=10)
    ax.grid(True,alpha=0.6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    ax.set_ylabel(f" tilt angle [\u00b0]", fontsize=20)
    ax.set_xlabel(f"time [UTC]", fontsize=20)
    ax.tick_params("x", rotation=45) 
    ax.legend(bbox_to_anchor=(0.1, 1.02),ncols=3,fontsize=8)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/tilt_vs_time_{string}_{device}.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def plot_tilt_time_icm_ids_height_subplots(string,device,icm_list,moni_files,plotFolder,geometry_list):
    fig = plt.figure(figsize=(8,5*4))
    gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    for i,i_icm in enumerate(icm_list[:]):
        ax = fig.add_subplot(gs[i,0])
        time_list, tilt_list = get_time_tilt(string,i_icm,moni_files,geometry_list)
        print(f"icm id {i_icm} has {len(time_list)} tilt measurements")
        time_list = pd.to_datetime(time_list, format='ISO8601')
        ax.plot(time_list,tilt_list,'-o',markersize=2,label=f"ICM {i_icm}",alpha=0.5)
        # print(f"time list {time_list}")
        ax.tick_params(axis='both',which='both', direction='in', labelsize=5)
        ax.grid(True,alpha=0.6)
        ax.set_ylim(0,2.6)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
        ax.set_ylabel(f" tilt angle [\u00b0]", fontsize=5)
        ax.set_xlabel(f"time [UTC]", fontsize=5)
        ax.tick_params("x", rotation=45) 
        if i < len(icm_list)-1:
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_xlabel(None) 
        ax.legend(fontsize=5)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/tilt_vs_time_{string}_{device}_height.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def plot_tilt_time_icm_ids_height_subplots_large(string,device,icm_list,moni_files,plotFolder,geometry_list,string_map_list):
    fig = plt.figure(figsize=(8,5*len(icm_list)))
    gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.5)
    for i,i_icm in enumerate(icm_list[:]):
        ax = fig.add_subplot(gs[i,0])
        time_list, tilt_list = get_time_tilt(string,i_icm,moni_files,geometry_list)
        print(f"icm id {i_icm} has {len(time_list)} tilt measurements")
        time_list = pd.to_datetime(time_list, format='ISO8601')
        upgrade_string,port = get_string_port_from_icm_id(i_icm,string_map_list)
        depth = get_depth_from_icm_id(i_icm,geometry_list)
        position,upgrade_string2 = get_string_position_from_icm_id(i_icm,geometry_list)
        # ax.text(0.05, 0.95, f"S{upgrade_string}:{port} OM {int(position):d} {depth:.0f} m", transform=ax.transAxes, ha='left', fontsize=5)
        ax.text(0.05, 0.92, f"S{upgrade_string} port {port} OM {int(position):d}", transform=ax.transAxes, ha='left', fontsize=14)
        ax.plot(time_list,tilt_list,'-o',markersize=2,label=f"{'tilt'}",alpha=0.5)
        # print(f"time list {time_list}")
        ax.tick_params(axis='both',which='both', direction='in', labelsize=14)
        ax.grid(True,alpha=0.6)
        ax.set_ylim(0,2.6)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
        ax.set_ylabel(f" tilt angle [\u00b0]", fontsize=14)
        ax.set_xlabel(f"time [UTC]", fontsize=12)
        ax.tick_params("x", rotation=45) 
        # if i < len(icm_list)-1:
        #     plt.setp(ax.get_xticklabels(), visible=False)
        #     ax.set_xlabel(None) 
        ax.legend(fontsize=14)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/tilt_vs_time_{string}_{device}_height_large.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def plot_tilt_vs_depth_icm_(string,device,icm_list,moni_files,plotFolder,geometry_list):
    fig = plt.figure(figsize=(8,5*4))
    gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    for i,i_icm in enumerate(icm_list[:]):
        ax = fig.add_subplot(gs[i,0])
        time_list, tilt_list = get_time_tilt(string,device,i_icm,moni_files,geometry_list)
        print(f"icm id {i_icm} has {len(time_list)} tilt measurements")
        time_list = pd.to_datetime(time_list, format='ISO8601')
        ax.plot(time_list,tilt_list,'-o',markersize=2,label=f"ICM {i_icm}",alpha=0.5)
        # print(f"time list {time_list}")
        ax.tick_params(axis='both',which='both', direction='in', labelsize=5)
        ax.grid(True,alpha=0.6)
        ax.set_ylim(0,2.6)
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
        ax.set_ylabel(f" tilt angle [\u00b0]", fontsize=5)
        ax.set_xlabel(f"time [UTC]", fontsize=5)
        ax.tick_params("x", rotation=45) 
        if i < len(icm_list)-1:
            plt.setp(ax.get_xticklabels(), visible=False)
            ax.set_xlabel(None) 
        ax.legend(fontsize=5)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/tilt_vs_time_{string}_{device}_height.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def main():
    from sensor_measurement import xDOMSensorMeasurement

    plotFolder = home+"/research_ua/icecube/upgrade/upgrade_comissioning/plots/"
    moni_data_folder = "/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/data/moni-data/"
    moni_files = sorted(glob.glob(moni_data_folder+"/*.json"))
    monidaq_device_list_file = "/Users/epaudel/research_ua/icecube/upgrade/timing_calibration/scripts/monidaq_device_list.json"
    string_map_dir = "/Users/epaudel/Library/CloudStorage/OneDrive-TheUniversityofAlabama/research_notes/notes/string_maps/"
    string_map_list = sorted(glob.glob(string_map_dir+"/*.json"))
    geometry_folder = "/Users/epaudel/research_ua/icecube/software/upgrade_commissioning_scripts/geometry" #upgrade commissioning scripts repo
    geometry_list = sorted(glob.glob(geometry_folder+"/string_*geometry.json"))
    tilt_plot("92","mDOM",moni_files,plotFolder)
    # tilt_plot("92","DEgg",moni_files,plotFolder)
    #get list of all devices on fieldhub 92
    icm_ids_92 = get_icm_id_list("92",["degg","mdom"],geometry_list)
    print(f"icm ids on string 92: {len(get_icm_id_list("92",["degg","mdom"],geometry_list))}")
    print(f"icm ids on string 92: {len(get_icm_id_list("92",["mdom"],geometry_list))}")
    print(f"icm ids on string 92: {len(get_icm_id_list("92",["degg"],geometry_list))}")
    print(f"icm ids on string 92: {len(get_icm_id_list_monidaq('92',['degg','mdom'],monidaq_device_list_file))}")
    icm_ids_92_monidaq = get_icm_id_list_monidaq("92",["degg","mdom"],monidaq_device_list_file)
    icm_ids_91_monidaq = get_icm_id_list_monidaq("91",["degg","mdom"],monidaq_device_list_file)
    icm_ids_90_monidaq = get_icm_id_list_monidaq("90",["degg","mdom"],monidaq_device_list_file)
    icm_ids_89_monidaq = get_icm_id_list_monidaq("89",["degg","mdom"],monidaq_device_list_file)
    icm_ids_88_monidaq = get_icm_id_list_monidaq("88",["degg","mdom"],monidaq_device_list_file)
    # icm_ids_92_monidaq = get_icm_id_list_monidaq("92",["degg","mdom"],monidaq_device_list_file)
    icm_ids_92_monidaq_2150_2300 = [icm_id for icm_id in icm_ids_92_monidaq if 2150 <= get_depth_from_icm_id(icm_id,geometry_list) <= 2300]
    depths_92 = [get_depth_from_icm_id(icm_id,geometry_list) for icm_id in icm_ids_92_monidaq]
    for i_icm,idepth in zip(icm_ids_92_monidaq,depths_92):
        print(f"icm id {i_icm} has depth {idepth:.0f} m")

    # print(get_icm_id_list("92",["degg","mdom"],geometry_list))
    # tilt_plot_port("92","mDOM",[5154],moni_files,plotFolder)
    tilt_plot_port_range("92","mDOM",5098,5175,moni_files,plotFolder)
    icm_ids_92 = get_icm_id_list("92",["mdom"],geometry_list)
    # plot_tilt_time_icm_ids("92","mDOM",icm_ids_92_monidaq,moni_files,plotFolder,geometry_list)
    # plot_tilt_time_icm_ids("92","mDOM",icm_ids_92_monidaq_2150_2300,moni_files,plotFolder,geometry_list)
    # plot_tilt_time_icm_ids_height_subplots("92","mDOM",icm_ids_92_monidaq_2150_2300,moni_files,plotFolder,geometry_list)
    plot_tilt_time_icm_ids_height_subplots_large("92","mDOM",icm_ids_92,moni_files,plotFolder,geometry_list,string_map_list)
    # plot_mean_tilt_depth("92",icm_ids_92_monidaq,moni_files,plotFolder,geometry_list,date_start="2026-03-28T00:00:00.000000+00:00",date_end="2026-04-28T00:00:00.000000+00:00")
    # plot_mean_tilt_depth("91",icm_ids_91_monidaq,moni_files,plotFolder,geometry_list,date_start="2026-03-28T00:00:00.000000+00:00",date_end="2026-04-28T00:00:00.000000+00:00")
    # plot_mean_tilt_depth("90",icm_ids_90_monidaq,moni_files,plotFolder,geometry_list,date_start="2026-03-28T00:00:00.000000+00:00",date_end="2026-04-28T00:00:00.000000+00:00")
    # plot_mean_tilt_depth("89",icm_ids_89_monidaq,moni_files,plotFolder,geometry_list,date_start="2026-03-28T00:00:00.000000+00:00",date_end="2026-04-28T00:00:00.000000+00:00")
    # plot_mean_tilt_depth("88",icm_ids_88_monidaq,moni_files,plotFolder,geometry_list,date_start="2026-03-28T00:00:00.000000+00:00",date_end="2026-04-28T00:00:00.000000+00:00")


if __name__ == "__main__":
    main()