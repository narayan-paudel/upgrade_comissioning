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


# Eliminated the air gap between fiber optic and IDQ sensor window Work on data of June 19 2025, calculating delay between falling edge of trigger and rising edge of photodiode pulse
data_folder = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/"
data_folder_flasher = data_folder+"flasher_setup/June232025/"
which_chain = 1
meas_files_00806 = sorted(glob.glob(data_folder_flasher+f"Flash{which_chain}FDC00806*"))
meas_files_00840 = sorted(glob.glob(data_folder_flasher+f"Flash{which_chain}FDC00840*"))
meas_files_00890 = sorted(glob.glob(data_folder_flasher+f"Flash{which_chain}FDC00890*"))
meas_files_00963 = sorted(glob.glob(data_folder_flasher+f"Flash{which_chain}FDC00963*"))

print(meas_files_00806)



import pandas as pd


def extractCSV(filepath):
    df = pd.read_csv(filepath,header=0,skiprows=[0])
    rotations = df.columns.values
    return np.asarray(df[['Start_to_Stop1']].values).T[0]

def removeOutliers(x_list,low_lim,high_lim):
    return np.asarray([ix for ix in x_list if ix <= high_lim and ix >= low_lim])


def plot_MSP430ReadingHist(which_chain):
    for FDC in [806,840,890,963][:]:
        fig = plt.figure(figsize=(8,5))
        gs = gridspec.GridSpec(nrows=1,ncols=1)
        ax = fig.add_subplot(gs[0])
        bins=np.linspace(17,28,45)
        print(bins)
        # bins=np.linspace(0,220000,2201)
        for n,flasher in enumerate([1,2,3,4,5]):
            flasher_data = extractCSV(data_folder_flasher+f"Flash{which_chain}FDC00{FDC}LED{flasher}.csv")
            ax.hist(flasher_data,bins=bins,histtype="step",linewidth=2.5,color=colorsCustom[n], label=f"LED {flasher}",alpha=1)
            # ax.hist(flasher_data,histtype="step",linewidth=2.5, label=f"FDC{FDC}_{flasher}",alpha=0.8)
        ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
        ax.set_ylabel(r"count", fontsize=22)
        # ax.set_xlabel(r"Start - Stop [ns]", fontsize=22)
        ax.set_xlabel(r"time delay [ns]", fontsize=22)
        ax.grid(True,alpha=0.6)
        ax.set_ylim(0,130)
        # ax.set_xlim(17,2)
        ax.set_xticks([17,19,21,23,25,27])
        # ax.set_ylim(285,310)
        # ax.set_yscale("log")
        ax.legend(title=f"Daisy chain {FDC}",title_fontsize=13,fontsize=11,ncols=5)
        # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
        plt.savefig(plotFolder+f"/Flash{which_chain}FDC{FDC}FlasherDelayHist.png",transparent=False,bbox_inches='tight')
        plt.savefig(plotFolder+f"/Flash{which_chain}FDC{FDC}FlasherDelayHist.pdf",transparent=False,bbox_inches='tight')
        plt.close()

# plot_MSP430ReadingHist(which_chain=1)

def plot_MSP430ReadingMeans(which_chain):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    for n,FDC in enumerate([806,840,890,963][:]):
        print(f"FDC {FDC}")
        mean_list = []
        std_err_list = []
        for flasher in [1,2,3,4,5]:
            print(f"flasher {flasher}")
            flasher_data = removeOutliers(extractCSV(data_folder_flasher+f"Flash{which_chain}FDC00{FDC}LED{flasher}.csv"),0,500)
            means = np.mean(flasher_data)
            mean_list.append(means)
            std_err = np.std(flasher_data)/np.sqrt(len(flasher_data)-1)
            std_err_list.append(std_err)
        ax.errorbar([1,2,3,4,5],mean_list,yerr=std_err_list,fmt="-o",color=colorsCustom[n],lw=3.5,ms=8.5,label=f"FDC {FDC}",alpha=1)
        print(f"{FDC} mean list {mean_list} diff {np.diff(mean_list)} range {max(mean_list)-min(mean_list):.2f} average {(max(mean_list)-min(mean_list))/4.0:.2f}")
    mean_list = []
    std_err_list = []
    for flasher in [1,2,3,4,5]:
        print(f"FDC {FDC}")
        flasher_data_cumulative = []
        for FDC in [806,840,890,963]:
            print(f"flasher {flasher}")
            flasher_data = removeOutliers(extractCSV(data_folder_flasher+f"Flash{which_chain}FDC00{FDC}LED{flasher}.csv"),0,500)
            flasher_data_cumulative += list(flasher_data)
        means = np.mean(flasher_data_cumulative)
        mean_list.append(means)
        std_err = np.std(flasher_data_cumulative)/np.sqrt(len(flasher_data_cumulative)-1)
        std_err_list.append(std_err)
    print(f"mean list {mean_list} diff {np.diff(mean_list)} range {max(mean_list)-min(mean_list):.2f} average {(max(mean_list)-min(mean_list))/4.0:.2f}")
    ax.errorbar([1,2,3,4,5],mean_list,yerr=std_err_list,fmt="-o",color="gray",lw=3.5,ms=8.5,label=f"combined",alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" flasher", fontsize=22)
    ax.set_ylabel(r"time delay [ns]", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.set_xlim(0,100)
    # ax.set_ylim(280,290)
    ax.legend(ncols=2,fontsize=16)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/Flasher{which_chain}DelayMeans.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/Flasher{which_chain}DelayMeans.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_MSP430ReadingMeans(which_chain=1)




#LED1 before and after breakout board Sep 15 2025
data_folder = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/"
data_folder_flasher = data_folder+"flasher_setup/Sep152025/"
which_chain = 1
meas_files_00890 = sorted(glob.glob(data_folder_flasher+f"Flash{which_chain}FDC00890*"))

print(meas_files_00890)


def plot_MSP430ReadingHistBreakoutSwap(which_chain):
    for FDC in [890][:]:
        fig = plt.figure(figsize=(8,5))
        gs = gridspec.GridSpec(nrows=1,ncols=1)
        ax = fig.add_subplot(gs[0])
        bins=np.linspace(12,19,29)
        print(bins)
        # bins=np.linspace(0,220000,2201)
        for n,flasher in enumerate([1]):
            settings = ["BreakoutAfterLED1","BreakoutAtStart"]
            for isettings,setting in enumerate(settings):
                flasher_data = removeOutliers(extractCSV(data_folder_flasher+f"Flash{which_chain}FDC00{FDC}LED{flasher}{setting}.csv"),0,500)
                print(f"FDC{FDC} LED{flasher} {setting} mean {np.mean(flasher_data)} std {np.std(flasher_data)}")
                ax.hist(flasher_data,bins=bins,histtype="step",linewidth=2.5,color=colorsCustom[n], label=f"LED {flasher} {setting} ({np.mean(flasher_data):.1f} $\pm$ {np.std(flasher_data):.1f} ns)",alpha=1)
                n += 1
                # ax.hist(flasher_data,histtype="step",linewidth=2.5, label=f"FDC{FDC}_{flasher}",alpha=0.8)
        ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
        ax.set_ylabel(r"count", fontsize=22)
        # ax.set_xlabel(r"Start - Stop [ns]", fontsize=22)
        ax.set_xlabel(r"time delay [ns]", fontsize=22)
        ax.grid(True,alpha=0.6)
        ax.set_ylim(0,450)
        # ax.set_xlim(17,2)
        ax.set_xticks([12,14,16,18,20])
        # ax.set_ylim(285,310)
        # ax.set_yscale("log")
        ax.legend(loc="upper right",title=f"Daisy chain {FDC}",title_fontsize=11,fontsize=11,ncols=1)
        # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
        plt.savefig(plotFolder+f"/Flash{which_chain}FDC{FDC}FlasherDelayHistBreakoutSwap.png",transparent=False,bbox_inches='tight')
        plt.savefig(plotFolder+f"/Flash{which_chain}FDC{FDC}FlasherDelayHistBreakoutSwap.pdf",transparent=False,bbox_inches='tight')
        plt.close()

plot_MSP430ReadingHistBreakoutSwap(which_chain=1)

def plot_MSP430ReadingMeans(which_chain):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    for n,FDC in enumerate([890][:]):
        print(f"FDC {FDC}")
        mean_list = []
        std_err_list = []
        for flasher in [1]:
            print(f"flasher {flasher}")
            flasher_data = removeOutliers(extractCSV(data_folder_flasher+f"Flash{which_chain}FDC00{FDC}LED{flasher}.csv"),0,500)
            means = np.mean(flasher_data)
            mean_list.append(means)
            std_err = np.std(flasher_data)/np.sqrt(len(flasher_data)-1)
            std_err_list.append(std_err)
        ax.errorbar([1,2,3,4,5],mean_list,yerr=std_err_list,fmt="-o",color=colorsCustom[n],lw=3.5,ms=8.5,label=f"FDC {FDC}",alpha=1)
        print(f"{FDC} mean list {mean_list} diff {np.diff(mean_list)} range {max(mean_list)-min(mean_list):.2f} average {(max(mean_list)-min(mean_list))/4.0:.2f}")
    mean_list = []
    std_err_list = []
    for flasher in [1,2,3,4,5]:
        print(f"FDC {FDC}")
        flasher_data_cumulative = []
        for FDC in [806,840,890,963]:
            print(f"flasher {flasher}")
            flasher_data = removeOutliers(extractCSV(data_folder_flasher+f"Flash{which_chain}FDC00{FDC}LED{flasher}.csv"),0,500)
            flasher_data_cumulative += list(flasher_data)
        means = np.mean(flasher_data_cumulative)
        mean_list.append(means)
        std_err = np.std(flasher_data_cumulative)/np.sqrt(len(flasher_data_cumulative)-1)
        std_err_list.append(std_err)
    print(f"mean list {mean_list} diff {np.diff(mean_list)} range {max(mean_list)-min(mean_list):.2f} average {(max(mean_list)-min(mean_list))/4.0:.2f}")
    ax.errorbar([1,2,3,4,5],mean_list,yerr=std_err_list,fmt="-o",color="gray",lw=3.5,ms=8.5,label=f"combined",alpha=1)
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" flasher", fontsize=22)
    ax.set_ylabel(r"time delay [ns]", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.set_xlim(0,100)
    # ax.set_ylim(280,290)
    ax.legend(ncols=2,fontsize=16)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/Flasher{which_chain}DelayMeans.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/Flasher{which_chain}DelayMeans.pdf",transparent=False,bbox_inches='tight')
    plt.close()

# plot_MSP430ReadingMeans(which_chain=1)
