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
data_folder_flasher = data_folder+"flasher_setup/June242025/"
which_chain = 1
# meas_files = {"No filter, DAC = 60000":data_folder_flasher+"/Flash1FDC00890LED1WOFilterDAC60000.csv", "With filter, DAC = 60000":data_folder_flasher+"/Flash1FDC00890LED1WFilterDAC60000.csv", "With filter, DAC = 55000":data_folder_flasher+"/Flash1FDC00890LED1WFilterDAC55000.csv"}
meas_files = {"With filter, DAC = 60000":data_folder_flasher+"/Flash1FDC00890LED1WFilterDAC60000.csv", "With filter, DAC = 55000":data_folder_flasher+"/Flash1FDC00890LED1WFilterDAC55000.csv"}

import pandas as pd


def extractCSV(filepath):
    df = pd.read_csv(filepath,header=0,skiprows=[0])
    rotations = df.columns.values
    return np.asarray(df[['Start_to_Stop1']].values).T[0]

def removeOutliers(x_list,low_lim,high_lim):
    return np.asarray([ix for ix in x_list if ix <= high_lim and ix >= low_lim])


def plot_MSP430ReadingHist(which_chain):
    for FDC in [890][:]:
        fig = plt.figure(figsize=(8,5))
        gs = gridspec.GridSpec(nrows=1,ncols=1)
        ax = fig.add_subplot(gs[0])
        bins=np.linspace(17,35,19)
        # print(bins)
        # bins=np.linspace(0,220000,2201)
        for n,isetting in enumerate(meas_files.keys()):
            flasher_data = extractCSV(meas_files[isetting])
            ax.hist(flasher_data,bins=bins,histtype="step",linewidth=2.5,color=colorsCustom[n], label=f"LED 1 {isetting}",alpha=1)
            # ax.hist(flasher_data,histtype="step",linewidth=2.5,color=colorsCustom[n], label=f"LED 1 {isetting}",alpha=1)
            # ax.hist(flasher_data,histtype="step",linewidth=2.5, label=f"FDC{FDC}_{flasher}",alpha=0.8)
        ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
        ax.set_ylabel(r"count", fontsize=22)
        # ax.set_xlabel(r"Start - Stop [ns]", fontsize=22)
        ax.set_xlabel(r"time delay [ns]", fontsize=22)
        ax.grid(True,alpha=0.6)
        # ax.set_ylim(0,130)
        # ax.set_xlim(17,2)
        # ax.set_xticks([17,19,21,23,25,27])
        # ax.set_ylim(285,310)
        # ax.set_yscale("log")
        ax.legend(title=f"Daisy chain {FDC}",title_fontsize=13,fontsize=11,ncols=1)
        # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
        plt.savefig(plotFolder+f"/Flash{which_chain}FDC{FDC}FlasherDelayHistFilter.png",transparent=False,bbox_inches='tight')
        plt.savefig(plotFolder+f"/Flash{which_chain}FDC{FDC}FlasherDelayHistFilter.pdf",transparent=False,bbox_inches='tight')
        plt.close()

plot_MSP430ReadingHist(which_chain=1)

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
