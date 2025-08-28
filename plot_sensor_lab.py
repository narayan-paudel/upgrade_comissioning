#!/usr/bin/env python3

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# sensor_file = "/home/enpaudel/icecube/upgrade/NTSFlashers/data/sensor_data/sensorsm184_20250725.csv"

path_to_sensor_files = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/"
dirs = ["North","East", "South","West"]
meas_day = "20250729"
mb = "mDOMMB"
file_list = [path_to_sensor_files+"sensors"+mb+"_"+idir+"_"+meas_day+".txt" for idir in dirs]
print(file_list)

plotFolder = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/plots/"


def read_sensor_file(sensor_file):
    with open(sensor_file) as txt_file:
        next(txt_file)
        gx = []
        gy = []
        gz = []
        bx = []
        by = []
        bz = []
        for row in txt_file:
            irow = row.split(" ")
            cols = [float(imeas) for imeas in irow]
            gx.append(cols[0])
            gy.append(cols[1])
            gz.append(cols[2])
            bx.append(cols[3])
            by.append(cols[4])
            bz.append(cols[5])
    return gx,gy,gz,bx,by,bz


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

def spherical_lists(x_list,y_list,z_list):
    r_list = []
    θ_list = []
    φ_list = []
    for x,y,z in zip(x_list,y_list,z_list):
        r,theta,phi = to_spherical(x,y,z)
        r_list.append(r)
        θ_list.append(theta)
        φ_list.append(phi)
    return r_list,θ_list,φ_list

# def get_mean_orientation(sensor_file):





def plot_ginclinations_mdoms():
    fig = plt.figure(figsize=(8,5))
    # print(time_stamps)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0])
    shift = 0
    xtick_positions = []
    xtick_labels = []
    for isensor, ipath in zip(dirs,file_list):
        gx,gy,gz,bx,by,bz = read_sensor_file(ipath)
        r,theta,phi = spherical_lists(gx,gy,gz)
        print(f"{isensor} {np.rad2deg(theta)}")
        x_list = [np.sin(itheta) for itheta in theta]
        y_list = [np.cos(itheta) for itheta in theta]

        mean_x = np.mean(x_list)
        mean_y = np.mean(y_list)
        mean_xy = np.sqrt(mean_x**2+mean_y**2)
        mean_x/=mean_xy
        mean_y/=mean_xy
        print(mean_x,mean_y)
        # ax.plot([0,mean_bx],[0,mean_by],"-",alpha=1)
        ax.arrow(0+shift,0,mean_x,mean_y,width=0.01,alpha=1)
        ax.plot(0+shift,0,"o",ms=8,alpha=1)
        xtick_positions.append(0+shift)
        xtick_labels.append(isensor)
        shift += 1
    
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
    plt.savefig(plotFolder+f"/LabSensorgInclinationmDOMs.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/LabSensorgInclinationmDOMs.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_ginclinations_mdoms()


def plot_binclinations_mdoms():
    fig = plt.figure(figsize=(8,5))
    # print(time_stamps)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0])
    shift = 0
    xtick_positions = []
    xtick_labels = []
    for isensor, ipath in zip(dirs,file_list):
        gx,gy,gz,bx,by,bz = read_sensor_file(ipath)
        r,theta,phi = spherical_lists(bx,by,bz)
        print(f"{isensor} {np.rad2deg(theta)}")
        x_list = [np.sin(itheta) for itheta in theta]
        y_list = [np.cos(itheta) for itheta in theta]

        mean_x = np.mean(x_list)
        mean_y = np.mean(y_list)
        mean_xy = np.sqrt(mean_x**2+mean_y**2)
        mean_x/=mean_xy
        mean_y/=mean_xy
        print(mean_x,mean_y)
        # ax.plot([0,mean_bx],[0,mean_by],"-",alpha=1)
        ax.arrow(0+shift,0,mean_x,mean_y,width=0.01,alpha=1)
        ax.plot(0+shift,0,"o",ms=8,alpha=1)
        xtick_positions.append(0+shift)
        xtick_labels.append(isensor)
        shift += 1
    
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
    plt.savefig(plotFolder+f"/LabSensorbInclinationmDOMs.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/LabSensorbInclinationmDOMs.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_binclinations_mdoms()

def plot_gorientations_mdoms():
    fig = plt.figure(figsize=(8,5))
    # print(time_stamps)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0])
    shift = 0
    xtick_positions = []
    xtick_labels = []
    for isensor, ipath in zip(dirs,file_list):
        gx,gy,gz,bx,by,bz = read_sensor_file(ipath)
        mean_bx = np.mean(gx)
        mean_by = np.mean(gy)
        mean_bxby = np.sqrt(mean_bx**2+mean_by**2)
        mean_bx/=mean_bxby
        mean_by/=mean_bxby
        mean_bz = np.mean(bz)
        print(mean_bx,mean_by)
        # ax.plot([0,mean_bx],[0,mean_by],"-",alpha=1)
        ax.arrow(0+shift,0,mean_bx,mean_by,width=0.01,alpha=1)
        ax.plot(0+shift,0,"o",ms=8,alpha=1)
        xtick_positions.append(0+shift)
        xtick_labels.append(isensor)
        shift += 1
    
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
    plt.savefig(plotFolder+f"/LabSensorgOrientationsmDOMs.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/LabSensorgOrientationsmDOMs.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_gorientations_mdoms()


def plot_borientations_mdoms():
    fig = plt.figure(figsize=(8,5))
    # print(time_stamps)
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0])
    shift = 0
    xtick_positions = []
    xtick_labels = []
    for isensor, ipath in zip(dirs,file_list):
        gx,gy,gz,bx,by,bz = read_sensor_file(ipath)
        mean_bx = np.mean(bx)
        mean_by = np.mean(by)
        mean_bxby = np.sqrt(mean_bx**2+mean_by**2)
        mean_bx/=mean_bxby
        mean_by/=mean_bxby
        mean_bz = np.mean(bz)
        print(mean_bx,mean_by)
        # ax.plot([0,mean_bx],[0,mean_by],"-",alpha=1)
        ax.arrow(0+shift,0,mean_bx,mean_by,width=0.01,alpha=1)
        ax.plot(0+shift,0,"o",ms=8,alpha=1)
        xtick_positions.append(0+shift)
        xtick_labels.append(isensor)
        shift += 1
    
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
    plt.savefig(plotFolder+f"/LabSensorbOrientationsmDOMs.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/LabSensorbOrientationsmDOMs.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_borientations_mdoms()