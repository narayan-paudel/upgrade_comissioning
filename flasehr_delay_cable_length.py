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

colorsCustom = ['#8dd3c7','#bebada','#80b1d3','#fdb462','#b3de69','#fb8072']



# Work on data of June 16 2025
data_folder = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/"
data_folder_flasher = data_folder+"flasher_setup/June182025/"


import pandas as pd


def extractCSV(filepath):
    # print(f"filepath {filepath}")
    df = pd.read_csv(filepath,header=0,skiprows=[0])
    rotations = df.columns.values
    return np.asarray(df[['Start_to_Stop1']].values).T[0]

def removeOutliers(x_list,low_lim,high_lim):
    return np.asarray([ix for ix in x_list if ix <= high_lim and ix >= low_lim])

FDC_list = [890]
LED_list = [1,2,3,4,5]
# coax_list = {1:38.2,2:156.0,3:188.0,4:344.0,5:1565.7,6:-119.3}
coax_list = {1:0.382,2:1.560,3:1.880,4:3.440,5:15.657}

def plot_MSP430ReadingHist():
    for FDC in FDC_list[:]:
        fig = plt.figure(figsize=(8,5))
        gs = gridspec.GridSpec(nrows=1,ncols=1)
        ax = fig.add_subplot(gs[0])
        # bins=np.linspace(285,310,51)
        bins=np.linspace(200,310,221)
        # bins=np.linspace(0,220000,2201)
        for n,cable in enumerate(coax_list.keys()):
            flasher_data = extractCSV(data_folder_flasher+f"/Flash1FDC00{FDC}LED1Coax{cable}.csv")
            ax.hist(flasher_data,bins=bins,histtype="step",linewidth=2.5,color=colorsCustom[n], label=f"cable length {coax_list[cable]:.3f} m",alpha=1)
            # ax.hist(flasher_data,histtype="step",linewidth=2.5, label=f"FDC{FDC}_{flasher}",alpha=0.8)
        ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
        ax.set_ylabel(r"count", fontsize=22)
        # ax.set_xlabel(r"Start - Stop [ns]", fontsize=22)
        ax.set_xlabel(r"time delay [ns]", fontsize=22)
        ax.grid(True,alpha=0.6)
        # ax.set_ylim(0,85)
        # ax.set_xlim(284,311)
        # ax.set_ylim(285,310)
        # ax.set_yscale("log")
        ax.legend()
        # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
        plt.savefig(plotFolder+f"/{FDC}_CableDelayhist.png",transparent=False,bbox_inches='tight')
        plt.savefig(plotFolder+f"/{FDC}_CableDelayhist.pdf",transparent=False,bbox_inches='tight')
        plt.close()

plot_MSP430ReadingHist()

def plot_MSP430ReadingMeans():
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    for n,FDC in enumerate(FDC_list[:]):
        print(f"FDC {FDC}")
        mean_list = []
        std_err_list = []
        for cable in coax_list.keys():
            print(f"flasher {cable}")
            flasher_data = removeOutliers(extractCSV(data_folder_flasher+f"/Flash1FDC00{FDC}LED1Coax{cable}.csv"),0,500)
            means = np.mean(flasher_data)
            mean_list.append(means)
            std_err = np.std(flasher_data)/np.sqrt(len(flasher_data)-1)
            std_err_list.append(std_err)
        ax.errorbar(coax_list.values(),mean_list,yerr=std_err_list,fmt="-o",color=colorsCustom[n],lw=3.5,ms=8.5,label=f"FDC_{FDC}",alpha=1)
        print(f"{FDC} mean list {mean_list} range {max(mean_list)-min(mean_list):.2f} distance range {max(coax_list.values())-min(coax_list.values())} rate speed {(max(coax_list.values())-min(coax_list.values()))/((max(mean_list)-min(mean_list))):.2f}")
    mean_list = []
    std_err_list = []
    # for cable in coax_list.keys():
    #     print(f"FDC {FDC}")
    #     flasher_data_cumulative = []
    #     for FDC in FDC_list:
    #         print(f"cable {cable}")
    #         flasher_data = removeOutliers(extractCSV(data_folder_flasher+f"/Flash1FDC00{FDC}LED1Coax{cable}.csv"),0,500)
    #         flasher_data_cumulative += list(flasher_data)
    #     means = np.mean(flasher_data_cumulative)
    #     mean_list.append(means)
    #     std_err = np.std(flasher_data_cumulative)/np.sqrt(len(flasher_data_cumulative)-1)
    #     std_err_list.append(std_err)
    # print(f"mean list {mean_list} diff {np.diff(mean_list)} range {max(mean_list)-min(mean_list):.2f} average {(max(mean_list)-min(mean_list))/4.0:.2f}")
    # ax.errorbar(coax_list.values(),mean_list,yerr=std_err_list,fmt="-o",color="gray",lw=3.5,ms=8.5,label=f"combined",alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" cable length [m]", fontsize=22)
    ax.set_ylabel(r"time delay [ns]", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.set_xlim(0,100)
    # ax.set_ylim(32,34)
    ax.legend(ncols=2,fontsize=18)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"cableDelayMeans_limits.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"cableDelayMeans_limits.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_MSP430ReadingMeans()

gap_status = {"AirGap":"with_air_gap","noAirGap":"no_air_gap"}
def plot_MSP430ReadingHistAirGap():
    for FDC in FDC_list[:]:
        fig = plt.figure(figsize=(8,5))
        gs = gridspec.GridSpec(nrows=1,ncols=1)
        ax = fig.add_subplot(gs[0])
        # bins=np.linspace(285,310,51)
        bins=np.linspace(285,310,51)
        # bins=np.linspace(0,220000,2201)
        for n,gap in enumerate(gap_status.keys()):
            flasher_data = extractCSV(data_folder_flasher+f"/Flash1FDC00{FDC}LED1{gap}.csv")
            ax.hist(flasher_data,bins=bins,histtype="step",linewidth=2.5,color=colorsCustom[n], label=f"{gap_status[gap]}",alpha=1)
            # ax.hist(flasher_data,histtype="step",linewidth=2.5, label=f"FDC{FDC}_{flasher}",alpha=0.8)
        ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
        ax.set_ylabel(r"count", fontsize=22)
        # ax.set_xlabel(r"Start - Stop [ns]", fontsize=22)
        ax.set_xlabel(r"time delay [ns]", fontsize=22)
        ax.grid(True,alpha=0.6)
        # ax.set_ylim(0,85)
        # ax.set_xlim(284,311)
        # ax.set_ylim(285,310)
        # ax.set_yscale("log")
        ax.legend()
        # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
        plt.savefig(plotFolder+f"/{FDC}_AirGaphist.png",transparent=False,bbox_inches='tight')
        plt.savefig(plotFolder+f"/{FDC}_AirGaphist.pdf",transparent=False,bbox_inches='tight')
        plt.close()

# plot_MSP430ReadingHistAirGap()


def plot_MSP430ReadingAirGapMeans():
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    for n,FDC in enumerate(FDC_list[:]):
        print(f"FDC {FDC}")
        mean_list = []
        std_err_list = []
        for gap in gap_status.keys():
            flasher_data = removeOutliers(extractCSV(data_folder_flasher+f"/Flash1FDC00{FDC}LED1{gap}.csv"),0,500)
            means = np.mean(flasher_data)
            mean_list.append(means)
            std_err = np.std(flasher_data)/np.sqrt(len(flasher_data)-1)
            std_err_list.append(std_err)
        ax.errorbar(gap_status.values(),mean_list,yerr=std_err_list,fmt="-o",color=colorsCustom[n],lw=3.5,ms=8.5,label=f"FDC_{FDC}",alpha=1)
        print(f"{FDC} mean list {mean_list} diff {np.diff(mean_list)} range {max(mean_list)-min(mean_list):.2f} average {(max(mean_list)-min(mean_list))/4.0:.2f}")
    mean_list = []
    std_err_list = []
    # for gap in gap_status.keys():
    #     print(f"FDC {FDC}")
    #     flasher_data_cumulative = []
    #     for FDC in FDC_list:
    #         # print(f"cable {cable}")
    #         flasher_data = removeOutliers(extractCSV(data_folder_flasher+f"/Flash1FDC00{FDC}LED1{gap}.csv"),0,500)
    #         flasher_data_cumulative += list(flasher_data)
    #     means = np.mean(flasher_data_cumulative)
    #     mean_list.append(means)
    #     std_err = np.std(flasher_data_cumulative)/np.sqrt(len(flasher_data_cumulative)-1)
    #     std_err_list.append(std_err)
    # print(f"mean list {mean_list} diff {np.diff(mean_list)} range {max(mean_list)-min(mean_list):.2f} average {(max(mean_list)-min(mean_list))/4.0:.2f}")
    # ax.errorbar(gap_status.values(),mean_list,yerr=std_err_list,fmt="-o",color="gray",lw=3.5,ms=8.5,label=f"combined",alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" flasher", fontsize=22)
    ax.set_ylabel(r"time delay [ns]", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.set_xlim(0,100)
    # ax.set_ylim(32,34)
    ax.legend(ncols=2,fontsize=18)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"airGapMeans_limits.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"airGapMeans_limits.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_MSP430ReadingAirGapMeans()