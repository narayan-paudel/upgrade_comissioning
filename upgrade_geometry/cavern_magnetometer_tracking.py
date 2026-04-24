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

from util import tilt_angle, to_spherical, to_spherical_lists


geo_to_IC0 = 1948.07 #distance from surface of IC geometry to IC z=0

from tilt_monitor_cavern import (get_icm_id_list, get_icm_id_list_monidaq,
get_depth_from_icm_id, get_mean_tilt_icm_id,
plot_tilt_time_icm_ids,plot_tilt_time_icm_ids_height_subplots, plot_mean_tilt_depth)


from geometry_util import (get_device_type, get_string_port_from_icm_id,
                            get_depth_from_icm_id, get_string_position_from_icm_id
                            ,get_icm_id_from_string_port)



def extract_float(line):
    pattern = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    numbers = re.findall(pattern, line)
    return [float(num) for num in numbers[-4:]]

def get_block(file_path, start_string, end_string):
    blocks = []
    current_block = []
    inside_block = False
    with open(file_path, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if start_string in line:
            inside_block = True
            current_block = []
            continue
        if inside_block:
            current_block.append(line)
        if end_string in line and inside_block:
            current_block.append(line)
            blocks.append(current_block)
            inside_block = False
            continue
    return blocks

def extract_data_from_blocks(blocks):
    data_list = []
    for block in blocks:
        data_block = {}
        for line in block:
            if "Start time" in line:
                # print(f"time line {line}")
                date_time_str = re.search(r'(\d{4}-\d{2}-\d{2}) \d{2}:\d{2}:(\d{2})', line)
                date = date_time_str.group(1)
                data_block["date"] = date
                continue
            if "Magnetometer/Accelerometer read successful" in line:
                string_port_str = re.search(r'fieldhub(\d{2}):(\d{4})', line)
                data_block["string"] = string_port_str.group(1)
                data_block["port"] = string_port_str.group(2)
            if " magnetometer_read: {" in line:
                bx,by,bz,t_mag = extract_float(line)
                data_block["bx"] = bx*10**6 #convert to microtesla
                data_block["by"] = by*10**6 #convert to microtesla
                data_block["bz"] = bz*10**6 #convert to microtesla
                data_block["t_mag"] = t_mag
            if "accelerometer_read: {" in line:
                gx,gy,gz,t_acc = extract_float(line)
                data_block["gx"] = gx
                data_block["gy"] = gy
                data_block["gz"] = gz
                data_block["t_acc"] = t_acc
        data_list.append(data_block)
        # print(f"Extracted data for block: {data_block}")
    return data_list


def extract_blocks(file_path, search_string, block_size):
    results = []

    with open(file_path, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if search_string in line:
            data_block = {}
            # print(f"start time {lines[i+1:i+2]}")
            block = lines[i+1:i+1+block_size]
            for iline in block:
                if "Start time" in iline:
                    date_time_str = re.search(r'(\d{4}-\d{2}-\d{2}) \d{2}:\d{2}:(\d{2})', iline)
                    date = date_time_str.group(1)
                    data_block["date"] = date
                    continue
            for iline in block:            
                if "Magnetometer/Accelerometer read successful" in iline:
                    string_port_str = re.search(r'fieldhub(\d{2}):(\d{4})', iline)
                    # print(f"string_port: {string_port_str.group(1)}, port: {string_port_str.group(2)}")
                    data_block["string"] = string_port_str.group(1)
                    data_block["port"] = string_port_str.group(2)


                # print(f"iline ({iline})")
                if " magnetometer_read: {" in iline:
                    # print(f"iline {iline}")
                    iline = str(iline)
                    bx,by,bz,t = extract_float(iline)
                    # print(f"bx {bx}, by {by}, bz {bz}")
                    data_block["bx"] = bx*10**6 #convert to microtesla
                    data_block["by"] = by*10**6 #convert to microtesla
                    data_block["bz"] = bz*10**6 #convert to microtesla
                    data_block["t_mag"] = t
                if "accelerometer_read: {" in iline:
                    # print(f"iline {iline}")
                    iline = str(iline)
                    gx,gy,gz,t = extract_float(iline)
                    # print(f"gx {gx}, gy {gy}, gz {gz}")
                    data_block["gx"] = gx
                    data_block["gy"] = gy
                    data_block["gz"] = gz
                    data_block["t_acc"] = t
            results.append(data_block)

    return results


def get_magnetometer_cartesian(measurement_list, string, port):
    # meas_this_port = [type(imeas) for imeas in measurement_list]
    # meas_this_port = [imeas for imeas in measurement_list]
    # print(f"meas this port {meas_this_port}")
    # for imeas in measurement_list:
    #     print(f"measurement {imeas}")
    #     print(f"imeas {imeas['string']} {imeas['port']}")
    meas_this_port = [imeas for imeas in measurement_list if imeas["string"] == string and imeas["port"] == port]
    print(f"meas this port {meas_this_port}")
    dates = [imeas["date"] for imeas in meas_this_port]
    date_list = list(set(dates))
    print(f"date list {date_list}")
    bx_mean_list = []
    bx_std_list = []
    by_mean_list = []
    by_std_list = []
    bz_mean_list = []
    bz_std_list = []
    for idate in date_list:
        bx_list = [imeas["bx"] for imeas in meas_this_port if imeas["date"] == idate]
        by_list = [imeas["by"] for imeas in meas_this_port if imeas["date"] == idate]
        bz_list = [imeas["bz"] for imeas in meas_this_port if imeas["date"] == idate]
        bx_mean = np.mean(bx_list)
        bx_std = np.std(bx_list)
        by_mean = np.mean(by_list)
        by_std = np.std(by_list)
        bz_mean = np.mean(bz_list)
        bz_std = np.std(bz_list)
        bx_mean_list.append(bx_mean)
        bx_std_list.append(bx_std)
        by_mean_list.append(by_mean)
        by_std_list.append(by_std)
        bz_mean_list.append(bz_mean)
        bz_std_list.append(bz_std)
    return date_list, bx_mean_list, bx_std_list, by_mean_list, by_std_list, bz_mean_list, bz_std_list


def get_temperature(measurement_list, string, port):
    meas_this_port = [imeas for imeas in measurement_list if imeas["string"] == string and imeas["port"] == port]
    print(f"meas this port {meas_this_port}")
    dates = [imeas["date"] for imeas in meas_this_port]
    date_list = list(set(dates))
    print(f"date list {date_list}")
    temp_B_mean_list  = []
    temp_B_std_list = []
    temp_g_mean_list  = []
    temp_g_std_list = []
    for idate in date_list:
        temp_B_list = [imeas["t_mag"] for imeas in meas_this_port if imeas["date"] == idate]
        temp_g_list = [imeas["t_acc"] for imeas in meas_this_port if imeas["date"] == idate]
        temp_B_mean_list.append(np.mean(temp_B_list))
        temp_B_std_list.append(np.std(temp_B_list))
        temp_g_mean_list.append(np.mean(temp_g_list))
        temp_g_std_list.append(np.std(temp_g_list))
    return date_list, temp_B_mean_list, temp_B_std_list, temp_g_mean_list, temp_g_std_list



def get_magnetometer_spherical(measurement_list, string, port):
    meas_this_port = [imeas for imeas in measurement_list if imeas["string"] == string and imeas["port"] == port]
    dates = [imeas["date"] for imeas in meas_this_port]
    date_list = list(set(dates))
    B_mean_list = []
    B_std_list = []
    theta_B_mean_list = []
    theta_B_std_list = []
    phi_B_mean_list = []
    phi_B_std_list = []
    for idate in date_list:
        bx_list = [imeas["bx"] for imeas in meas_this_port if imeas["date"] == idate]
        by_list = [imeas["by"] for imeas in meas_this_port if imeas["date"] == idate]
        bz_list = [imeas["bz"] for imeas in meas_this_port if imeas["date"] == idate]
        B, theta_B, phi_B = to_spherical_lists(bx_list, by_list, bz_list)
        B_mean = np.mean(B)
        B_std = np.std(B)
        theta_B_mean = np.mean(theta_B)
        theta_B_std = np.std(theta_B)
        phi_B_mean = np.mean(phi_B)
        phi_B_std = np.std(phi_B)
        B_mean_list.append(B_mean)
        B_std_list.append(B_std)
        theta_B_mean_list.append(theta_B_mean)
        theta_B_std_list.append(theta_B_std)
        phi_B_mean_list.append(phi_B_mean)
        phi_B_std_list.append(phi_B_std)
    return date_list, B_mean_list, B_std_list, theta_B_mean_list, theta_B_std_list, phi_B_mean_list, phi_B_std_list

def plot_Bx_By_Bz(time_list,Bx,Bx_std,By,By_std,Bz,Bz_std,plotFolder,plot_suffix,select_string, select_port,B_ylim=None):
    fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig,hspace=0.5)
    ax = fig.add_subplot(gs[0,0])
    time_list = pd.to_datetime(time_list, format='ISO8601')
    print(f"len(time_list), len(Bx), len(By), len(Bz): {time_list}, {Bx}, {By}, {Bz}")
    #sorting by time
    sorted_lists = zip(*sorted(zip(time_list, Bx, Bx_std,By, By_std, Bz, Bz_std)))
    time_list, Bx, Bx_std, By, By_std, Bz, Bz_std = map(list, sorted_lists)
    ax.text(0.05, 0.95, f"S{select_string} port {select_port}", transform=ax.transAxes, ha='left', fontsize=14)
    ax.errorbar(time_list, Bx, yerr=Bx_std, fmt='-o', markersize=4, label=f"{r'B$_{{x}}$'}", alpha=1)
    ax.errorbar(time_list, By, yerr=By_std, fmt='-o', markersize=4, label=f"{r'B$_{{y}}$'}", alpha=1)
    ax.errorbar(time_list, Bz, yerr=Bz_std, fmt='-o', markersize=4, label=f"{r'B$_{{z}}$'}", alpha=1)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=14)
    ax.grid(True,alpha=0.6)
    # ax.set_ylim(0,2.6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    ax.set_ylabel(r"B$_{i}$ [$\mu$T]", fontsize=14)
    ax.set_xlabel(f"time [UTC]", fontsize=12)
    ax.tick_params("x", rotation=45) 
    # if i < len(icm_list)-1:
    #     plt.setp(ax.get_xticklabels(), visible=False)
    #     ax.set_xlabel(None) 
    ax.legend(fontsize=14)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/BxByBz_vs_time_{plot_suffix}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_B(time_list,B,B_std,plotFolder,plot_suffix,select_string, select_port, B_ylim=None):
    fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig,hspace=0.5)
    ax = fig.add_subplot(gs[0,0])
    time_list = pd.to_datetime(time_list, format='ISO8601')
    print(f"len(time_list), len(B), len(B_std): {time_list}, {B}, {B_std}")
    #sorting by time
    sorted_lists = zip(*sorted(zip(time_list, B, B_std)))
    time_list, B, B_std = map(list, sorted_lists)
    ax.text(0.05, 0.95, f"S{select_string} port {select_port}", transform=ax.transAxes, ha='left', fontsize=14)
    ax.errorbar(time_list, B, yerr=B_std, fmt='-o', markersize=4, label=f"{r'B'}", alpha=1)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=14)
    ax.grid(True,alpha=0.6)
    # ax.set_ylim(0,2.6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    ax.set_ylabel(r"B [$\mu$T]", fontsize=14)
    ax.set_xlabel(f"time [UTC]", fontsize=14)
    ax.tick_params("x", rotation=45) 
    if B_ylim is not None:
        ax.set_ylim(*B_ylim)
    # if i < len(icm_list)-1:
    #     plt.setp(ax.get_xticklabels(), visible=False)
    #     ax.set_xlabel(None) 
    ax.legend(fontsize=14)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/B_vs_time_{plot_suffix}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_B_theta(time_list,theta_B_mean,theta_B_std,plotFolder,plot_suffix,select_string, select_port,B_theta_ylim=None):
    fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig,hspace=0.5)
    ax = fig.add_subplot(gs[0,0])
    time_list = pd.to_datetime(time_list, format='ISO8601')
    theta_B_mean = [np.rad2deg(itheta) for itheta in theta_B_mean]
    theta_B_std = [np.rad2deg(istd) for istd in theta_B_std]
    #sorting by time
    sorted_lists = zip(*sorted(zip(time_list, theta_B_mean, theta_B_std)))
    time_list, theta_B_mean, theta_B_std = map(list, sorted_lists)
    ax.text(0.05, 0.95, f"S{select_string} port {select_port}", transform=ax.transAxes, ha='left', fontsize=14)
    ax.errorbar(time_list, theta_B_mean, yerr=theta_B_std, fmt='-o', markersize=4, label=f"{r'$\theta$ [°]'}", alpha=1)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=14)
    ax.grid(True,alpha=0.6)
    # ax.set_ylim(0,2.6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    ax.set_ylabel(r"$\theta$ [°]", fontsize=14)
    ax.set_xlabel(f"time [UTC]", fontsize=14)
    ax.tick_params("x", rotation=45) 
    # if i < len(icm_list)-1:
    #     plt.setp(ax.get_xticklabels(), visible=False)
    #     ax.set_xlabel(None) 
    ax.legend(fontsize=14)
    if B_theta_ylim is not None:
        ax.set_ylim(*B_theta_ylim)
    # ax.set_ylim(0,20)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/B_theta_vs_time_{plot_suffix}.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def plot_B_phi(time_list,phi_B_mean,phi_B_std,plotFolder,plot_suffix,select_string, select_port,B_phi_ylim=None):
    fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig,hspace=0.5)
    ax = fig.add_subplot(gs[0,0])
    time_list = pd.to_datetime(time_list, format='ISO8601')
    phi_B_mean = [np.rad2deg(iphi) for iphi in phi_B_mean]
    phi_B_std = [np.rad2deg(istd) for istd in phi_B_std]
    #sorting by time
    sorted_lists = zip(*sorted(zip(time_list, phi_B_mean, phi_B_std)))
    time_list, phi_B_mean, phi_B_std = map(list, sorted_lists)
    ax.text(0.05, 0.95, f"S{select_string} port {select_port}", transform=ax.transAxes, ha='left', fontsize=14)
    ax.errorbar(time_list, phi_B_mean, yerr=phi_B_std, fmt='-o', markersize=4, label=f"{r'$\phi$ [°]'}", alpha=1)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=14)
    ax.grid(True,alpha=0.6)
    # ax.set_ylim(0,2.6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    ax.set_ylabel(r"$\phi$ [°]", fontsize=14)
    ax.set_xlabel(f"time [UTC]", fontsize=14)
    if B_phi_ylim is not None:
        ax.set_ylim(*B_phi_ylim)
    ax.tick_params("x", rotation=45) 
    # if i < len(icm_list)-1:
    #     plt.setp(ax.get_xticklabels(), visible=False)
    #     ax.set_xlabel(None) 
    ax.legend(fontsize=14)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/B_phi_vs_time_{plot_suffix}.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def get_accelerometer_cartesian(measurement_list, string, port):
    meas_this_port = [imeas for imeas in measurement_list if imeas["string"] == string and imeas["port"] == port]
    print(f"meas this port {meas_this_port}")
    dates = [imeas["date"] for imeas in meas_this_port]
    date_list = list(set(dates))
    print(f"date list {date_list}")
    gx_mean_list = []
    gx_std_list = []
    gy_mean_list = []
    gy_std_list = []
    gz_mean_list = []
    gz_std_list = []
    for idate in date_list:
        gx_list = [imeas["gx"] for imeas in meas_this_port if imeas["date"] == idate]
        gy_list = [imeas["gy"] for imeas in meas_this_port if imeas["date"] == idate]
        gz_list = [imeas["gz"] for imeas in meas_this_port if imeas["date"] == idate]
        gx_mean = np.mean(gx_list)
        gx_std = np.std(gx_list)
        gy_mean = np.mean(gy_list)
        gy_std = np.std(gy_list)
        gz_mean = np.mean(gz_list)
        gz_std = np.std(gz_list)
        gx_mean_list.append(gx_mean)
        gx_std_list.append(gx_std)
        gy_mean_list.append(gy_mean)
        gy_std_list.append(gy_std)
        gz_mean_list.append(gz_mean)
        gz_std_list.append(gz_std)
    return date_list, gx_mean_list, gx_std_list, gy_mean_list, gy_std_list, gz_mean_list, gz_std_list




def get_accelerometer_spherical(measurement_list, string, port):
    meas_this_port = [imeas for imeas in measurement_list if imeas["string"] == string and imeas["port"] == port]
    dates = [imeas["date"] for imeas in meas_this_port]
    date_list = list(set(dates))
    g_mean_list = []
    g_std_list = []
    theta_g_mean_list = []
    theta_g_std_list = []
    phi_g_mean_list = []
    phi_g_std_list = []
    for idate in date_list:
        gx_list = [imeas["gx"] for imeas in meas_this_port if imeas["date"] == idate]
        gy_list = [imeas["gy"] for imeas in meas_this_port if imeas["date"] == idate]
        gz_list = [imeas["gz"] for imeas in meas_this_port if imeas["date"] == idate]
        g, theta_g, phi_g = to_spherical_lists(gx_list, gy_list, gz_list)
        g_mean = np.mean(g)
        g_std = np.std(g)
        theta_g_mean = np.mean(theta_g)
        theta_g_std = np.std(theta_g)
        phi_g_mean = np.mean(phi_g)
        phi_g_std = np.std(phi_g)
        g_mean_list.append(g_mean)
        g_std_list.append(g_std)
        theta_g_mean_list.append(theta_g_mean)
        theta_g_std_list.append(theta_g_std)
        phi_g_mean_list.append(phi_g_mean)
        phi_g_std_list.append(phi_g_std)
    return date_list, g_mean_list, g_std_list, theta_g_mean_list, theta_g_std_list, phi_g_mean_list, phi_g_std_list

def plot_gx_gy_gz(time_list,gx,gx_std,gy,gy_std,gz,gz_std,plotFolder,plot_suffix,select_string, select_port,g_ylim=None):
    fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig,hspace=0.5)
    ax = fig.add_subplot(gs[0,0])
    time_list = pd.to_datetime(time_list, format='ISO8601')
    print(f"len(time_list), len(gx), len(gy), len(gz): {time_list}, {gx}, {gy}, {gz}")
    #sorting by time
    sorted_lists = zip(*sorted(zip(time_list, gx, gy, gz, gx_std, gy_std, gz_std)))
    time_list, gx, gy, gz, gx_std, gy_std, gz_std = map(list, sorted_lists)
    ax.text(0.05, 0.95, f"S{select_string} port {select_port}", transform=ax.transAxes, ha='left', fontsize=14)
    ax.errorbar(time_list, gx, yerr=gx_std, fmt='-o', markersize=4, label=f"{r'g$_{{x}}$'}", alpha=1)
    ax.errorbar(time_list, gy, yerr=gy_std, fmt='-o', markersize=4, label=f"{r'g$_{{y}}$'}", alpha=1)
    ax.errorbar(time_list, gz, yerr=gz_std, fmt='-o', markersize=4, label=f"{r'g$_{{z}}$'}", alpha=1)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=14)
    ax.grid(True,alpha=0.6)
    # ax.set_ylim(0,2.6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    ax.set_ylabel(r"g$_{i}$ [ms$^{-2}$]", fontsize=14)
    ax.set_xlabel(f"time [UTC]", fontsize=12)
    ax.tick_params("x", rotation=45) 
    # if i < len(icm_list)-1:
    #     plt.setp(ax.get_xticklabels(), visible=False)
    #     ax.set_xlabel(None) 
    ax.legend(fontsize=14)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/gxgygz_vs_time_{plot_suffix}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_g(time_list,g,g_std,plotFolder,plot_suffix,select_string, select_port,g_ylim=None):
    fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig,hspace=0.5)
    ax = fig.add_subplot(gs[0,0])
    time_list = pd.to_datetime(time_list, format='ISO8601')
    #sorting by time
    sorted_lists = zip(*sorted(zip(time_list, g, g_std)))
    time_list, g, g_std = map(list, sorted_lists)
    ax.text(0.05, 0.95, f"S{select_string} port {select_port}", transform=ax.transAxes, ha='left', fontsize=14)
    ax.errorbar(time_list, g, yerr=g_std, fmt='-o', markersize=4, label=f"{r'g'}", alpha=1)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=14)
    ax.grid(True,alpha=0.6)
    # ax.set_ylim(0,2.6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    ax.set_ylabel(r"g [ms$^{-2}$]", fontsize=14)
    ax.set_xlabel(f"time [UTC]", fontsize=14)
    ax.tick_params("x", rotation=45) 
    if g_ylim is not None:
        ax.set_ylim(*g_ylim)
    # if i < len(icm_list)-1:
    #     plt.setp(ax.get_xticklabels(), visible=False)
    #     ax.set_xlabel(None) 
    ax.legend(fontsize=14)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/g_vs_time_{plot_suffix}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_temperature(time_list, temp_B_mean_list, temp_B_std_list, temp_g_mean_list, temp_g_std_list,plotFolder,plot_suffix,select_string, select_port, temp_ylim=None):
    fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig,hspace=0.5)
    ax = fig.add_subplot(gs[0,0])
    time_list = pd.to_datetime(time_list, format='ISO8601')
    #sorting by time
    sorted_lists = zip(*sorted(zip(time_list, temp_B_mean_list, temp_B_std_list, temp_g_mean_list, temp_g_std_list)))
    time_list, temp_B_mean_list, temp_B_std_list, temp_g_mean_list, temp_g_std_list = map(list, sorted_lists)
    ax.text(0.05, 0.95, f"S{select_string} port {select_port}", transform=ax.transAxes, ha='left', fontsize=14)
    ax.errorbar(time_list, temp_B_mean_list, yerr=temp_B_std_list, fmt='-o', markersize=4, label=f"{r'Temp B'}", alpha=1)
    ax.errorbar(time_list, temp_g_mean_list, yerr=temp_g_std_list, fmt='-s', markersize=4, label=f"{r'Temp g'}", alpha=1)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=14)
    ax.grid(True,alpha=0.6)
    # ax.set_ylim(0,2.6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    ax.set_ylabel(r"Temperature [°C]", fontsize=14)
    ax.set_xlabel(f"time [UTC]", fontsize=14)
    ax.tick_params("x", rotation=45) 
    if temp_ylim is not None:
        ax.set_ylim(*temp_ylim)
    # if i < len(icm_list)-1:
    #     plt.setp(ax.get_xticklabels(), visible=False)
    #     ax.set_xlabel(None) 
    ax.legend(fontsize=14)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/temperature_vs_time_{plot_suffix}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_g_theta(time_list,theta_g_mean,theta_g_std,plotFolder,plot_suffix,select_string, select_port,g_theta_ylim=None):
    fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig,hspace=0.5)
    ax = fig.add_subplot(gs[0,0])
    time_list = pd.to_datetime(time_list, format='ISO8601')
    theta_g_mean = [np.rad2deg(itheta) for itheta in theta_g_mean]
    theta_g_std = [np.rad2deg(istd) for istd in theta_g_std]
    #sorting by time
    sorted_lists = zip(*sorted(zip(time_list, theta_g_mean, theta_g_std)))
    time_list, theta_g_mean, theta_g_std = map(list, sorted_lists)
    ax.text(0.05, 0.95, f"S{select_string} port {select_port}", transform=ax.transAxes, ha='left', fontsize=14)
    ax.errorbar(time_list, theta_g_mean, yerr=theta_g_std, fmt='-o', markersize=4, label=f"{r'$\theta$ [°]'}", alpha=1)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=14)
    ax.grid(True,alpha=0.6)
    # ax.set_ylim(0,2.6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    ax.set_ylabel(r"$\theta$ [°]", fontsize=14)
    ax.set_xlabel(f"time [UTC]", fontsize=14)
    ax.tick_params("x", rotation=45) 
    # if i < len(icm_list)-1:
    #     plt.setp(ax.get_xticklabels(), visible=False)
    #     ax.set_xlabel(None) 
    ax.legend(fontsize=14)
    if g_theta_ylim is not None:
        ax.set_ylim(*g_theta_ylim)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/g_theta_vs_time_{plot_suffix}.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def plot_g_phi(time_list,phi_g_mean,phi_g_std,plotFolder,plot_suffix,select_string, select_port,g_phi_ylim=None):
    fig = plt.figure(figsize=(8,5))
    # gs = gridspec.GridSpec(nrows=len(icm_list),ncols=1, figure=fig,hspace=0.1)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig,hspace=0.5)
    ax = fig.add_subplot(gs[0,0])
    time_list = pd.to_datetime(time_list, format='ISO8601')
    phi_g_mean = [np.rad2deg(iphi) for iphi in phi_g_mean]
    phi_g_std = [np.rad2deg(istd) for istd in phi_g_std]
    #sorting by time
    sorted_lists = zip(*sorted(zip(time_list, phi_g_mean, phi_g_std)))
    time_list, phi_g_mean, phi_g_std = map(list, sorted_lists)
    ax.text(0.05, 0.95, f"S{select_string} port {select_port}", transform=ax.transAxes, ha='left', fontsize=14)
    ax.errorbar(time_list, phi_g_mean, yerr=phi_g_std, fmt='-o', markersize=4, label=f"{r'$\phi$ [°]'}", alpha=1)
    # print(f"time list {time_list}")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=14)
    ax.grid(True,alpha=0.6)
    # ax.set_ylim(0,2.6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d', tz=None))
    ax.set_ylabel(r"$\phi$ [°]", fontsize=14)
    ax.set_xlabel(f"time [UTC]", fontsize=14)
    if g_phi_ylim is not None:
        ax.set_ylim(*g_phi_ylim)
    ax.tick_params("x", rotation=45) 
    # if i < len(icm_list)-1:
    #     plt.setp(ax.get_xticklabels(), visible=False)
    #     ax.set_xlabel(None) 
    ax.legend(fontsize=14)
    # plt.savefig(plotFolder+f"/tilt_vs_time_{board_name}_{hostname}_{port}_{icm_id}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/g_phi_vs_time_{plot_suffix}.pdf",transparent=False,bbox_inches='tight')
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

    log_files_path = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/string_92_logs/"
    # log_files = glob.glob(log_files_path+"Camera-LOG-String92_FreezeInRuns_SpecialRun_*.txt")
    log_files = sorted(glob.glob(log_files_path+"04*/Camera-LOG-*.txt"))

    print(f"Number of log files found: {len(log_files)}")

    print(log_files)
    # Example usage
    measurement_list = []
    for ifile in log_files[:]:
        # print(f"Extracting blocks from {ifile}")
        # blocks = get_block(ifile, "magnetometer/accelerometer reading STARTS!!!", "All tasks completed successfully!")
        blocks = get_block(ifile, "magnetometer/accelerometer reading STARTS!!!", "Magnetometer/Accelerometer read TRIAL (host_port)")
        # print(blocks)
        print(f"Number of blocks extracted from {ifile}: {len(blocks)}")
        # blocks = extract_blocks(ifile, "magnetometer/accelerometer reading STARTS!!!", 10)
        data_list = extract_data_from_blocks(blocks)
        # print(f"data list from {ifile}: {data_list}")
        measurement_list.extend(data_list)
    # print(f"Total measurements collected: {measurement_list}")
    
    select_string_list = ["92"]
    # select_port_list = ["5108","5106"]
    select_port_list = ["5108","5106","5098","5095","5090"]
    y_adjustment = "standard"
    # y_adjustment = "none"
    if y_adjustment == "standard":
        g_ylim = {"5108":[0,15], "5106":[0,15], "5098":[0,15], "5095":[0,15], "5090":[0,15]}
        g_theta_ylim = {"5108": [0,2], "5106":[0,2], "5098":[0,2], "5095":[0,2], "5090":[0,2]}
        g_phi_ylim = {"5108": [140,200], "5106":[300,360], "5098":[50,110], "5095":[40,100], "5090":[300,360]}
        B_ylim = {"5108":[0,90], "5106":[0,90], "5098":[0,90], "5095":[0,90], "5090":[0,90]}
        B_theta_ylim = {"5108":[35,45], "5106":[0,2], "5098":[0,360], "5095":[25,35], "5090":[40,50]}
        B_phi_ylim = {"5108":[140,200], "5106":[300,360], "5098":[0,360], "5095":[300,360], "5090":[300,360]}
    elif y_adjustment == "full":
        g_ylim = {"5108":[0,15], "5106":[0,15], "5098":[0,15], "5095":[0,15], "5090":[0,15]}
        g_theta_ylim = {"5108": [0,90], "5106":[0,90], "5098":[0,90], "5095":[0,90], "5090":[0,90]}
        g_phi_ylim = {"5108": [0,360], "5106":[0,360], "5098":[0,360], "5095":[0,360], "5090":[0,360]}
        B_ylim = {"5108":[0,90], "5106":[0,90], "5098":[0,90], "5095":[0,90], "5090":[0,90]}
        B_theta_ylim = {"5108":[0,90], "5106":[0,90], "5098":[0,90], "5095":[0,90], "5090":[0,90]}
        B_phi_ylim = {"5108":[0,360], "5106":[0,360], "5098":[0,360], "5095":[0,360], "5090":[0,360]}
    elif y_adjustment == "none":
        g_ylim = {"5108":None, "5106":None, "5098":None, "5095":None, "5090":None}
        g_theta_ylim = {"5108": None, "5106":None, "5098":None, "5095":None, "5090":None}
        g_phi_ylim = {"5108": None, "5106":None, "5098":None, "5095":None, "5090":None}
        B_ylim = {"5108":None, "5106":None, "5098":None, "5095":None, "5090":None}
        B_theta_ylim = {"5108":None, "5106":None, "5098":None, "5095":None, "5090":None}
        B_phi_ylim = {"5108":None, "5106":None, "5098":None, "5095":None, "5090":None}

    # select_port_list = ["5106","5108","5300","5095","51００","5127"]
    # icm_id = get_icm_id_from_string_port(int(select_string), int(select_port), string_map_list)
    # string,position = get_string_position_from_icm_id(icm_id, geometry_list)

    for iport in select_port_list:
        for istring in select_string_list:
            select_string = istring
            select_port = iport
            date_list, Bx_mean, Bx_std, By_mean, By_std, Bz_mean, Bz_std = get_magnetometer_cartesian(measurement_list, select_string, select_port)
            date_list, B_mean, B_std, theta_B_mean, theta_B_std, phi_B_mean, phi_B_std = get_magnetometer_spherical(measurement_list, select_string, select_port)
            plot_Bx_By_Bz(date_list,Bx_mean,Bx_std,By_mean,By_std,Bz_mean,Bz_std,plotFolder,f"S{select_string}_p{select_port}_B_components",select_string, select_port, B_ylim=B_ylim[iport])
            plot_B(date_list,B_mean,B_std, plotFolder, f"S{select_string}_p{select_port}_B_components", select_string, select_port, B_ylim=B_ylim[iport])
            plot_B_theta(date_list,theta_B_mean,theta_B_std, plotFolder, f"S{select_string}_p{select_port}_B_components", select_string, select_port,B_theta_ylim=B_theta_ylim[iport])
            plot_B_phi(date_list,phi_B_mean,phi_B_std, plotFolder, f"S{select_string}_p{select_port}_B_components", select_string, select_port,B_phi_ylim=B_phi_ylim[iport])
            date_list, gx_mean, gx_std, gy_mean, gy_std, gz_mean, gz_std = get_accelerometer_cartesian(measurement_list, select_string, select_port)
            date_list, g_mean, g_std, theta_g_mean, theta_g_std, phi_g_mean, phi_g_std = get_accelerometer_spherical(measurement_list, select_string, select_port)
            plot_gx_gy_gz(date_list,gx_mean,gx_std,gy_mean,gy_std,gz_mean,gz_std,plotFolder,f"S{select_string}_p{select_port}_g_components",select_string, select_port,g_ylim=g_ylim[iport])
            plot_g(date_list,g_mean,g_std, plotFolder, f"S{select_string}_p{select_port}_g_components", select_string, select_port,g_ylim=g_ylim[iport])
            plot_g_theta(date_list,theta_g_mean,theta_g_std, plotFolder, f"S{select_string}_p{select_port}_g_components", select_string, select_port,g_theta_ylim=g_theta_ylim[iport])
            plot_g_phi(date_list,phi_g_mean,phi_g_std, plotFolder, f"S{select_string}_p{select_port}_g_components", select_string, select_port,g_phi_ylim=g_phi_ylim[iport])
            date_list, temp_B_mean_list, temp_B_std_list, temp_g_mean_list, temp_g_std_list = get_temperature(measurement_list, select_string, select_port)
            plot_temperature(date_list, temp_B_mean_list, temp_B_std_list, temp_g_mean_list, temp_g_std_list, plotFolder, f"S{select_string}_p{select_port}_temperature", select_string, select_port)
                    # plot_g(date_list,g_mean,g_std, plotFolder, f"S{select_string}_p{select_port}_g_components_zoomed", select_string, select_port,g_ylim=g_ylim[iport])
if __name__ == "__main__":
    main()
