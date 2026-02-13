#!/usr/bin/env python3

import glob
import json
import re

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pathlib import Path
home = str(Path.home())

from utils import to_spherical, tilt_angle, get_means, get_sub_df, to_spherical_list, to_spherical_list_rotated,meas_from_df
from calibrate_magnetometer_2d import corrected_ellipse

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']


plotFolder = home+"/research_ua/icecube/upgrade/upgrade_comissioning/plots/"
# print(degg_list)

string_maps_folder = "/Users/epaudel/Library/CloudStorage/OneDrive-TheUniversityofAlabama/research_notes/notes/string_maps/"
string_maps_list = sorted(glob.glob(string_maps_folder+"*.json"))
string_maps_map = {}
for istring_map in string_maps_list:
    string_number = int(re.findall(r'-?\d*\.?\d+', istring_map.split("/")[-1].split("_")[1])[0])
    string_maps_map[string_number] = istring_map
print(string_maps_map)
ICU_sensor_measurements = "/Users/epaudel/Library/CloudStorage/OneDrive-TheUniversityofAlabama/research_notes/notes/sensor_measurements_at_ICU/"

geometry_cmd = "/Users/epaudel/Library/CloudStorage/OneDrive-TheUniversityofAlabama/research_notes/notes/geometry/Upgrade_Strings_CMD-v103_final_final_final.xlsx"

meas_file_list = sorted(glob.glob(ICU_sensor_measurements+"*.txt"))

meas_file = meas_file_list[-1]

print(f"meas file {meas_file}")

import pandas as pd

df = pd.read_csv(meas_file,header=0,sep=" ")
print(df.head())

def get_mean_g(df):
    '''
    returns field in ms^-2 and degrees
    '''
    gx,gy,gz = df[["gx","gy","gz"]].values.T
    r,theta,phi = to_spherical_list(gx,gy,gz)
    return np.mean(gx),np.std(gx),np.mean(gy),np.std(gy),np.mean(gz),np.std(gz),np.mean(r),np.std(r),np.rad2deg(np.mean(theta)),\
        np.rad2deg(np.std(theta)),np.rad2deg(np.mean(phi)),np.rad2deg(np.std(phi))


gx,gx_std,gy,gy_std,gz,gz_std,r,r_std,theta,theta_std,phi,phi_std = get_mean_g(df)

print(f"gx {gx:.2f} +- {gx_std:.2f} ms^-2 gy {gy:.2f} +- {gy_std:.2f} ms^-2 gz {gz:.2f} +- {gz_std:.2f} ms^-2 r {r:.2f} +- {r_std:.2f} ms^-2 theta {theta:.2f} +- {theta_std:.2f} deg phi {phi:.2f} +- {phi_std:.2f} deg")

gx_list = []
gy_list = []
gz_list = []
r_list = []
theta_list = []
phi_list = []

def plot_hist_g(x):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    bins = np.linspace(8.8,10.1,40)
    ax.hist(x,bins=bins,histtype="step",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.legend(,ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    # ax.set_ylim(-180,180)
    ax.set_ylabel(r" # of samples", fontsize=20)
    ax.set_xlabel(f"g [ms$^{-2}$]", fontsize=20)
    plt.savefig(plotFolder+f"/g_dist.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/g_dist.pdf",transparent=False,bbox_inches='tight')
    plt.close()

def plot_hist_phi(x):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # bins = np.linspace(8.8,10.1,40)
    # ax.hist(x,bins=bins,histtype="step",lw=2.5,alpha=1)
    ax.hist(x,histtype="step",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.legend(,ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    # ax.set_ylim(-180,180)
    ax.set_ylabel(r" # of samples", fontsize=20)
    ax.set_xlabel(f" $\\phi$\u00b0", fontsize=20)
    plt.savefig(plotFolder+f"/phi_dist.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/phi_dist.pdf",transparent=False,bbox_inches='tight')
    plt.close()


def plot_hist_theta(x):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # bins = np.linspace(-1,90,92)
    bins = np.linspace(-1.5,6.5,81)
    ax.hist(x,bins=bins,histtype="step",lw=2.5,alpha=1)
    # ax.hist(x,histtype="step",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.legend(,ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    # ax.set_ylim(-180,180)
    ax.set_ylabel(r" # of samples", fontsize=20)
    ax.set_xlabel(f" $\\theta$\u00b0", fontsize=20)
    plt.savefig(plotFolder+f"/theta_dist.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/theta_dist.pdf",transparent=False,bbox_inches='tight')
    plt.close()
    

for ielt in meas_file_list:
    df = pd.read_csv(ielt,header=0,sep=" ")
    gx,gx_std,gy,gy_std,gz,gz_std,r,r_std,theta,theta_std,phi,phi_std = get_mean_g(df)
    gx_list.append(gx)
    gy_list.append(gy)
    gz_list.append(gz)
    r_list.append(r)
    theta_list.append(theta)
    phi_list.append(phi)
    # print(f"file {ielt} gx {gx:.2f} +- {gx_std:.2f} ms^-2 gy {gy:.2f} +- {gy_std:.2f} ms^-2 gz {gz:.2f} +- {gz_std:.2f} ms^-2 r {r:.2f} +- {r_std:.2f} ms^-2 theta {theta:.2f} +- {theta_std:.2f} deg phi {phi:.2f} +- {phi_std:.2f} deg")
plot_hist_g(r_list)
plot_hist_phi(phi_list)
plot_hist_theta(theta_list)
theta_list = [i for i in theta_list if i>=-1 and i<=20]
print(f"mean theta {np.median(theta_list):.2f} deg std theta {np.std(theta_list):.2f} deg mean phi {np.mean(phi_list):.2f} deg std phi {np.std(phi_list):.2f} deg")


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


class xDOMSensorMeasurement:
    def __init__(self,df,string,port):
        self.df = df
        self.string = string
        self.port = port
        self.get_mean_B()
        self.get_mean_g()

    def get_mean_B(self):
        '''
        returns field in muT and degrees
        '''
        bx,by,bz = self.df[["bx","by","bz"]].values.T
        B,B_theta,B_phi = to_spherical_list(bx,by,bz)
        if np.rad2deg(np.mean(B_phi))<0:
            print(f"negative  phi detected {np.rad2deg(np.mean(B_phi))}")
        for ielt,iname in zip([bx,by,bz,B],["bx","by","bz","B"]):
            setattr(self, f"{iname}_mean", np.mean(ielt)*10**6)
            setattr(self, f"{iname}_std", np.std(ielt)*10**6)
        for ielt,iname in zip([B_theta,B_phi],["B_theta","B_phi"]):
            setattr(self, f"{iname}_mean", np.rad2deg(np.mean(ielt)))
            setattr(self, f"{iname}_std", np.rad2deg(np.std(ielt)))


    def get_mean_g(self):
        '''
        returns field in ms^-2 and degrees
        '''
        gx,gy,gz = self.df[["gx","gy","gz"]].values.T
        g,g_theta,g_phi = to_spherical_list(gx,gy,gz)
        if np.rad2deg(np.mean(g_phi))<0:
            print(f"negative  phi detected {np.rad2deg(np.mean(g_phi))}")
        for ielt,iname in zip([gx,gy,gz,g],["gx","gy","gz","g"]):
            setattr(self, f"{iname}_mean", np.mean(ielt))
            setattr(self, f"{iname}_std", np.std(ielt))
        for ielt,iname in zip([g_theta,g_phi],["g_theta","g_phi"]):
            setattr(self, f"{iname}_mean", np.rad2deg(np.mean(ielt)))
            setattr(self, f"{iname}_std", np.rad2deg(np.std(ielt)))

def xDOMSensorMeasurementList(meas_file_list):
    sensor_measurement_list = []
    for ifile in meas_file_list:
        string = int(re.findall(r'-?\d*\.?\d+', ifile.split("/")[-1].split("_")[1])[0])
        port = int(re.findall(r'-?\d*\.?\d+', ifile.split("/")[-1].split("_")[2])[0])
        # print(f"file {ifile} string {string} port {port}")
        df = pd.read_csv(ifile,header=0,sep=" ")
        sensor_measurement_list.append(xDOMSensorMeasurement(df,string,port))
    return sensor_measurement_list

sensor_measurement_list = xDOMSensorMeasurementList(meas_file_list)
print(f"no of sensor measurements {len(sensor_measurement_list)}")
for istring in string_maps_map.keys():
    sensor_measurement_istring = [i for i in sensor_measurement_list if i.string==istring]
    string_map = string_maps_map[istring]
    with open(string_map, 'r') as f:
        data = json.load(f)
    for isensor in sensor_measurement_istring:
        for imap in data:
            if imap["port"]==isensor.port:
                isensor.control_port = imap["control_port"]
                isensor.hostname = imap["hostname"]
                isensor.wp_addr = imap["wp_addr"]
                isensor.board_type_name = imap["board_type_name"]
                isensor.icm_id = imap["icm_id"]

def plot_tilt_angles(sensor_measurement_list):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    tilt_angles = []
    for isensor in sensor_measurement_list:
        tilt_angles.append(isensor.g_theta_mean)
    bins = np.linspace(0,90,91)
    ax.hist(tilt_angles,bins=bins,histtype="step",lw=2.5,alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    ax.set_ylabel(r" # of samples", fontsize=20)
    ax.set_xlabel(f"tilt angle [deg]", fontsize=20)
    plt.savefig(plotFolder+f"/tilt_angle_dist.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/tilt_angle_dist.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_tilt_angles(sensor_measurement_list)

for isensor in sensor_measurement_list:
    if isensor.g_theta_mean>70:
    print(f"{isensor.board_type_name}")