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

#measurement made on June 5 2025 for Daisy chain mDOM_FDC_00890
# data_folder = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/"
# flasher1 = data_folder+"flasher_setup/June62025/TextFileGraphflasher12.csv"
# flasher2 = data_folder+"flasher_setup/June62025/TextFileGraphflasher22.csv"
# flasher3 = data_folder+"flasher_setup/June62025/TextFileGraphflasher32.csv"
# flasher4 = data_folder+"flasher_setup/June62025/TextFileGraphflasher42.csv"
# flasher5 = data_folder+"flasher_setup/June62025/TextFileGraphflasher52.csv"

# #measurement made on June 10 2025 for Daisy chains mDOM_FDC_00890, MDOM_FDC_00806, MDOM_FDC_00840, MDOM_FDC_00963 
# #eg /Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/flasher_setup/June102025/TextFileGraph00806_flasher1.csv
# data_folder = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/"
# meas_files_00806 = sorted(glob.glob(data_folder+"flasher_setup/June102025/TextFileGraph00806*"))
# meas_files_00840 = sorted(glob.glob(data_folder+"flasher_setup/June102025/TextFileGraph00840*"))
# meas_files_00890 = sorted(glob.glob(data_folder+"flasher_setup/June102025/TextFileGraph00890*"))
# meas_files_00963 = sorted(glob.glob(data_folder+"flasher_setup/June102025/TextFileGraph00963*"))


# print(meas_files_00806)



# import pandas as pd


# def extractCSV(filepath):
#     df = pd.read_csv(filepath,header=0,skiprows=[0])
#     rotations = df.columns.values
#     return np.asarray(df[['Start_to_Stop1']].values).T[0]

# def removeOutliers(x_list,low_lim,high_lim):
#     return np.asarray([ix for ix in x_list if ix <= high_lim and ix >= low_lim])


# for FDC in [806,840,890,963]:
#     for flasher in [1,2,3,4,5]:
#         flasher_data = extractCSV(data_folder+f"flasher_setup/June102025/TextFileGraph00{FDC}_flasher{flasher}.csv")

# def plot_MSP430Reading(low_lim,high_lim):
#     for FDC in [806,840,890,963][:]:
#         fig = plt.figure(figsize=(8,5))
#         gs = gridspec.GridSpec(nrows=1,ncols=1)
#         ax = fig.add_subplot(gs[0])
#         for flasher in [1,2,3,4,5]:
#             flasher_data = removeOutliers(extractCSV(data_folder+f"flasher_setup/June102025/TextFileGraph00{FDC}_flasher{flasher}.csv"),low_lim,high_lim)
#             ax.plot(range(len(flasher_data)),flasher_data,"o",label=f"FDC{FDC}_{flasher}",alpha=0.8)
#         ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#         ax.set_xlabel(r" samples", fontsize=22)
#         ax.set_ylabel(r"Start - Stop [ns]", fontsize=22)
#         ax.grid(True,alpha=0.6)
#         ax.set_ylim(30,40)
#         # ax.set_ylim(30,40)
#         # ax.set_ylim(285,310)
#         ax.legend()
#         # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#         plt.savefig(plotFolder+f"/MSPReadingsJune10LightTightWFabricFDC{FDC}_limits{low_lim}_{high_lim}.png",transparent=False,bbox_inches='tight')
#         plt.savefig(plotFolder+f"/MSPReadingsJune10LightTightWFabricFDC{FDC}_limits{low_lim}_{high_lim}.pdf",transparent=False,bbox_inches='tight')
#         plt.close()

# # plot_MSP430Reading(30,40)

# def plot_MSP430ReadingMeans(low_lim,high_lim):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     for FDC in [806,840,890,963][:]:
#         print(f"FDC {FDC}")
#         mean_list = []
#         std_err_list = []
#         for flasher in [1,2,3,4,5]:
#             print(f"flasher {flasher}")
#             flasher_data = removeOutliers(extractCSV(data_folder+f"flasher_setup/June102025/TextFileGraph00{FDC}_flasher{flasher}.csv"),low_lim,high_lim)
#             means = np.mean(flasher_data)
#             mean_list.append(means)
#             std_err = np.std(flasher_data)/np.sqrt(len(flasher_data)-1)
#             std_err_list.append(std_err)
#         ax.errorbar([1,2,3,4,5],mean_list,yerr=std_err_list,fmt="-o",lw=3.5,ms=8.5,label=f"FDC_{FDC}",alpha=0.8)
#         print(f"{FDC} mean list {mean_list} diff {np.diff(mean_list)} range {max(mean_list)-min(mean_list):.2f} average {(max(mean_list)-min(mean_list))/4.0:.2f}")
#     mean_list = []
#     std_err_list = []
#     for flasher in [1,2,3,4,5]:
#         print(f"FDC {FDC}")
#         flasher_data_cumulative = []
#         for FDC in [806,840,890,963]:
#             print(f"flasher {flasher}")
#             flasher_data = removeOutliers(extractCSV(data_folder+f"flasher_setup/June102025/TextFileGraph00{FDC}_flasher{flasher}.csv"),low_lim,high_lim)
#             flasher_data_cumulative += list(flasher_data)
#         means = np.mean(flasher_data_cumulative)
#         mean_list.append(means)
#         std_err = np.std(flasher_data_cumulative)/np.sqrt(len(flasher_data_cumulative)-1)
#         std_err_list.append(std_err)
#     print(f"mean list {mean_list} diff {np.diff(mean_list)} range {max(mean_list)-min(mean_list):.2f} average {(max(mean_list)-min(mean_list))/4.0:.2f}")
#     ax.errorbar([1,2,3,4,5],mean_list,yerr=std_err_list,fmt="-o",color="gray",lw=3.5,ms=8.5,label=f"combined",alpha=0.8)
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#     ax.set_xlabel(r" flasher", fontsize=22)
#     ax.set_ylabel(r"delay with trigger [ns]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     # ax.set_xlim(0,100)
#     ax.set_ylim(32,34)
#     ax.legend(ncols=2,fontsize=18)
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/MSPReadingsJune10LightTightWFabricMeans_limits{low_lim}_{high_lim}.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/MSPReadingsJune10LightTightWFabricMeans_limits{low_lim}_{high_lim}.pdf",transparent=False,bbox_inches='tight')
#     plt.close()

# plot_MSP430ReadingMeans(30,40)












# flasher1_data = removeOutliers(extractCSV(flasher1))
# flasher2_data = removeOutliers(extractCSV(flasher2))
# flasher3_data = removeOutliers(extractCSV(flasher3))
# flasher4_data = removeOutliers(extractCSV(flasher4))
# flasher5_data = removeOutliers(extractCSV(flasher5))



# def plot_MSP430Reading(flasher1_data,flasher2_data,flasher3_data,flasher4_data,flasher5_data):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     ax.plot(range(len(flasher1_data)),flasher1_data,"-o",label="flasher_1")
#     ax.plot(range(len(flasher2_data)),flasher2_data,"-o",label="flasher_2")
#     ax.plot(range(len(flasher3_data)),flasher3_data,"-o",label="flasher_3")
#     ax.plot(range(len(flasher4_data)),flasher4_data,"-o",label="flasher_4")
#     ax.plot(range(len(flasher5_data)),flasher5_data,"-o",label="flasher_5")
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#     ax.set_xlabel(r" samples", fontsize=22)
#     ax.set_ylabel(r"Start - Stop [ns]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     # ax.set_xlim(0,100)
#     ax.set_ylim(285,310)
#     ax.legend()
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/MSPReadingsJune6LightTightWFabric.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/MSPReadingsJune6LightTightWFabric.pdf",transparent=False,bbox_inches='tight')
#     plt.close()

# plot_MSP430Reading(flasher1_data,flasher2_data,flasher3_data,flasher4_data,flasher5_data)


# def plot_MSP430ReadingMeans(flasher1_data,flasher2_data,flasher3_data,flasher4_data,flasher5_data):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     flasherList = [flasher1_data,flasher2_data,flasher3_data,flasher4_data,flasher5_data]
#     means = [np.mean(iflasher) for iflasher in flasherList]
#     print(f"max {max(means)} min {min(means)} diff {max(means)-min(means)} average {(max(means)-min(means))/4}")
#     stds = [np.std(iflasher) for iflasher in flasherList]
#     flasherLabels = [1,2,3,4,5]
#     ax.errorbar(flasherLabels,means,yerr=stds,fmt="-o",ms=8.5,label="flasher_1")
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#     ax.set_xlabel(r" flasher", fontsize=22)
#     ax.set_ylabel(r"delay with trigger [ns]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     # ax.set_xlim(0,100)
#     # ax.set_ylim(285,310)
#     # ax.legend()
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/MSPReadingsJune6LightTightWFabricMeans.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/MSPReadingsJune6LightTightWFabricMeans.pdf",transparent=False,bbox_inches='tight')
#     plt.close()

# plot_MSP430ReadingMeans(flasher1_data,flasher2_data,flasher3_data,flasher4_data,flasher5_data)









# def plot_MSP430FlasherReading(flasher_off,flasher_on):
#     fig = plt.figure(figsize=(8,5))
#     gs = gridspec.GridSpec(nrows=1,ncols=1)
#     ax = fig.add_subplot(gs[0])
#     ax.plot(range(len(flasher_off)),flasher_off,"-o",label="flasher off")
#     ax.plot(range(len(flasher_on)),flasher_on,"-o",label="flasher on")
#     ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
#     ax.set_xlabel(r" time", fontsize=22)
#     ax.set_ylabel(r"Start - Stop [ns]", fontsize=22)
#     ax.grid(True,alpha=0.6)
#     # ax.set_xlim(0,100)
#     ax.legend()
#     # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
#     plt.savefig(plotFolder+f"/FlasherReadings.png",transparent=False,bbox_inches='tight')
#     plt.savefig(plotFolder+f"/FlasherReadings.pdf",transparent=False,bbox_inches='tight')
#     plt.close()

# # IDQLightOffFlasherOff_list = extractCSV(IDQLightOffFlasherOff)
# # IDQLightOffFlasherON_list = extractCSV(IDQLightOffFlasherON)
# # plot_MSP430FlasherReading(IDQLightOffFlasherOff_list,IDQLightOffFlasherON_list)



# Work on data of June 16 2025
# data_folder = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/"
# meas_files_00806 = sorted(glob.glob(data_folder+"flasher_setup/June162025/Flash2FDC00806*"))
# meas_files_00840 = sorted(glob.glob(data_folder+"flasher_setup/June162025/Flash2FDC00840*"))
# meas_files_00890 = sorted(glob.glob(data_folder+"flasher_setup/June162025/Flash2FDC00890*"))
# meas_files_00963 = sorted(glob.glob(data_folder+"flasher_setup/June162025/Flash2FDC00963*"))

# # Eliminated the air gap between fiber optic and IDQ sensor window Work on data of June 19 2025
# data_folder = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/"
# data_folder_flasher = data_folder+"flasher_setup/June192025/"
# which_chain = 1
# meas_files_00806 = sorted(glob.glob(data_folder_flasher+f"Flash{which_chain}FDC00806*"))
# meas_files_00840 = sorted(glob.glob(data_folder_flasher+f"Flash{which_chain}FDC00840*"))
# meas_files_00890 = sorted(glob.glob(data_folder_flasher+f"Flash{which_chain}FDC00890*"))
# meas_files_00963 = sorted(glob.glob(data_folder_flasher+f"Flash{which_chain}FDC00963*"))

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
        bins=np.linspace(276,292,65)
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
        ax.set_ylim(0,81)
        ax.set_xlim(276,292)
        ax.set_xticks(np.linspace(278,292,8))
        # ax.set_ylim(285,310)
        # ax.set_yscale("log")
        ax.legend(title=f"Daisy chain {FDC}",title_fontsize=13,fontsize=13)
        # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
        plt.savefig(plotFolder+f"/Flash{which_chain}FDC{FDC}FlasherDelayHist.png",transparent=False,bbox_inches='tight')
        plt.savefig(plotFolder+f"/Flash{which_chain}FDC{FDC}FlasherDelayHist.pdf",transparent=False,bbox_inches='tight')
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
    ax.set_ylim(280,290)
    ax.legend(ncols=2,fontsize=16)
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/Flasher{which_chain}DelayMeans.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/Flasher{which_chain}DelayMeans.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_MSP430ReadingMeans(which_chain=1)
