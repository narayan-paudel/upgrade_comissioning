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

colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']


data_folder = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/"
siggen = data_folder+"flasher_setup/TextFileGraphWFGen.csv"
IDQLightOff = data_folder+"flasher_setup/TextFileGraphWithIDQLightOff.csv"
IDQLightOn = data_folder+"flasher_setup/TextFileGraphWithIDQLightOn.csv"
IDQLightOffFlasherOff = data_folder+"/flasher_setup/TextFileGraphWithADQLightOffFlasherOff.csv"
IDQLightOffFlasherON = data_folder+"/flasher_setup/TextFileGraphWithADQLightOffFlasherON.csv"

import pandas as pd


df = pd.read_csv(siggen,header=0,skiprows=[0])
print(df.columns.values)
TDC_values = np.asarray(df[['Start_to_Stop1']].values).T
# print(TDC_values)

def extractCSV(filepath):
    df = pd.read_csv(filepath,header=0,skiprows=[0])
    rotations = df.columns.values
    return np.asarray(df[['Start_to_Stop1']].values).T[0]
siggen_list = extractCSV(siggen)
print(siggen_list)
IDQLightOff_list = extractCSV(IDQLightOff)
IDQLightOn_list = extractCSV(IDQLightOn)

def plot_MSP430Reading(siggen_list,IDQLightOff_list,IDQLightOn_list):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.plot(range(len(siggen_list)),siggen_list,"-o",label="signal generator")
    ax.plot(range(len(IDQLightOff_list)),IDQLightOff_list,"-o",label="IDQ Light off")
    ax.plot(range(len(IDQLightOn_list)),IDQLightOn_list,"-o",label="IDQ Light on")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" time", fontsize=22)
    ax.set_ylabel(r"Start - Stop [ns]", fontsize=22)
    ax.grid(True,alpha=0.6)
    ax.set_xlim(0,100)
    ax.legend()
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/MSPReadings.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/MSPReadings.pdf",transparent=False,bbox_inches='tight')
    plt.close()

plot_MSP430Reading(siggen_list,IDQLightOff_list,IDQLightOn_list)

def plot_MSP430FlasherReading(flasher_off,flasher_on):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1)
    ax = fig.add_subplot(gs[0])
    ax.plot(range(len(flasher_off)),flasher_off,"-o",label="flasher off")
    ax.plot(range(len(flasher_on)),flasher_on,"-o",label="flasher on")
    ax.tick_params(axis='both',which='both', direction='in', labelsize=22)
    ax.set_xlabel(r" time", fontsize=22)
    ax.set_ylabel(r"Start - Stop [ns]", fontsize=22)
    ax.grid(True,alpha=0.6)
    # ax.set_xlim(0,100)
    ax.legend()
    # ax.legend(["B$_{x}$","B$_{y}$","B$_{z}$","B"],fontsize=14,ncols=4,bbox_to_anchor=(0.90, 0.95),loc="right")
    plt.savefig(plotFolder+f"/FlasherReadings.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/FlasherReadings.pdf",transparent=False,bbox_inches='tight')
    plt.close()

IDQLightOffFlasherOff_list = extractCSV(IDQLightOffFlasherOff)
IDQLightOffFlasherON_list = extractCSV(IDQLightOffFlasherON)
plot_MSP430FlasherReading(IDQLightOffFlasherOff_list,IDQLightOffFlasherON_list)