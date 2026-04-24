#!/usr/bin/env python

import os
import glob

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pathlib import Path
home = str(Path.home())


import pandas as pd
import json

def main():
    df_mDOM_pole_run1 = pd.read_csv(
        "/Users/epaudel/research_ua/icecube/upgrade/pole_rotation/sensorsmDOMMB_20251215_run1.txt",
        header=0,sep=" ")
    df_mDOM_pole_run2 = pd.read_csv(
        "/Users/epaudel/research_ua/icecube/upgrade/pole_rotation/sensorsmDOMMB_20251215_run2.txt",
        header=0,sep=" ")
    df_mDOM_pole_run3 = pd.read_csv(
        "/Users/epaudel/research_ua/icecube/upgrade/pole_rotation/sensorsmDOMMB_20251215_run3.txt",
        header=0,sep=" ")
    df_mDOM_pole_run4 = pd.read_csv(
        "/Users/epaudel/research_ua/icecube/upgrade/pole_rotation/sensorsmDOMMB_20251215_run4.txt",
        header=0,sep=" ")
    run = 4
    if run == 4:
        df_mDOM = df_mDOM_pole_run4
    elif run == 3:
        df_mDOM = df_mDOM_pole_run3    

    bx = df_mDOM["bx"].values
    by = df_mDOM["by"].values
    bz = df_mDOM["bz"].values
    gx = df_mDOM["gx"].values
    gy = df_mDOM["gy"].values
    gz = df_mDOM["gz"].values

    sensor_map = {}

    if run == 4:
        sensor_map["bx_PMT8"] = list(bx[200:225]) #PMT 8 North
        sensor_map["by_PMT8"] = list(by[200:225]) #PMT 8 North
        sensor_map["bz_PMT8"] = list(bz[200:225]) #PMT 8 North
        sensor_map["bx_PMT8p"] = list(bx[300:325]) #PMT 8 and 9 North
        sensor_map["by_PMT8p"] = list(by[300:325]) #PMT 8 and 9 North
        sensor_map["bz_PMT8p"] = list(bz[300:325]) 
        sensor_map["bx_PMT9"] = list(bx[400:425]) 
        sensor_map["by_PMT9"] = list(by[400:425]) 
        sensor_map["bz_PMT9"] = list(bz[400:425]) 
        sensor_map["bx_PMT9p"] = list(bx[500:525]) 
        sensor_map["by_PMT9p"] = list(by[500:525]) 
        sensor_map["bz_PMT9p"] = list(bz[500:525]) 
        sensor_map["bx_PMT10"] = list(bx[575:600]) 
        sensor_map["by_PMT10"] = list(by[575:600]) 
        sensor_map["bz_PMT10"] = list(bz[575:600]) 
        sensor_map["bx_PMT10p"] = list(bx[625:650]) 
        sensor_map["by_PMT10p"] = list(by[625:650]) 
        sensor_map["bz_PMT10p"] = list(bz[625:650]) 
        sensor_map["bx_PMT11"] = list(bx[700:725]) 
        sensor_map["by_PMT11"] = list(by[700:725]) 
        sensor_map["bz_PMT11"] = list(bz[700:725]) 
        sensor_map["bx_PMT11p"] = list(bx[750:775]) 
        sensor_map["by_PMT11p"] = list(by[750:775]) 
        sensor_map["bz_PMT11p"] = list(bz[750:775]) 
        sensor_map["bx_PMT4"] = list(bx[800:825]) 
        sensor_map["by_PMT4"] = list(by[800:825]) 
        sensor_map["bz_PMT4"] = list(bz[800:825]) 
        sensor_map["bx_PMT4p"] = list(bx[875:900]) 
        sensor_map["by_PMT4p"] = list(by[875:900]) 
        sensor_map["bz_PMT4p"] = list(bz[875:900]) 
        sensor_map["bx_PMT5"] = list(bx[925:950]) 
        sensor_map["by_PMT5"] = list(by[925:950]) 
        sensor_map["bz_PMT5"] = list(bz[925:950]) 
        sensor_map["bx_PMT5p"] = list(bx[975:1000]) 
        sensor_map["by_PMT5p"] = list(by[975:1000]) 
        sensor_map["bz_PMT5p"] = list(bz[975:1000]) 
        sensor_map["bx_PMT6"] = list(bx[1050:1075]) 
        sensor_map["by_PMT6"] = list(by[1050:1075]) 
        sensor_map["bz_PMT6"] = list(bz[1050:1075]) 
        sensor_map["bx_PMT6p"] = list(bx[1100:1125]) 
        sensor_map["by_PMT6p"] = list(by[1100:1125]) 
        sensor_map["bz_PMT6p"] = list(bz[1100:1125]) 
        sensor_map["bx_PMT7"] = list(bx[1150:1175]) 
        sensor_map["by_PMT7"] = list(by[1150:1175]) 
        sensor_map["bz_PMT7"] = list(bz[1150:1175]) 
        sensor_map["bx_PMT7p"] = list(bx[1225:1250]) 
        sensor_map["by_PMT7p"] = list(by[1225:1250]) 
        sensor_map["bz_PMT7p"] = list(bz[1225:1250]) 
        sensor_map["bx_PMT8"] = list(bx[1280:1305]) 
        sensor_map["by_PMT8"] = list(by[1280:1305]) 
        sensor_map["bz_PMT8"] = list(bz[1280:1305]) 
        #gravitational acceleration ma)p
        sensor_map["gx_PMT8"] = list(gx[200:225]) #PMT 8 Nort)h
        sensor_map["gy_PMT8"] = list(gy[200:225]) #PMT 8 Nort)h
        sensor_map["gz_PMT8"] = list(gz[200:225]) #PMT 8 Nort)h
        sensor_map["gx_PMT8p"] = list(gx[300:325]) #PMT 8 and 9 Nort)h
        sensor_map["gy_PMT8p"] = list(gy[300:325]) #PMT 8 and 9 Nort)h
        sensor_map["gz_PMT8p"] = list(gz[300:325]) 
        sensor_map["gx_PMT9"] = list(gx[400:425]) 
        sensor_map["gy_PMT9"] = list(gy[400:425]) 
        sensor_map["gz_PMT9"] = list(gz[400:425]) 
        sensor_map["gx_PMT9p"] = list(gx[500:525]) 
        sensor_map["gy_PMT9p"] = list(gy[500:525]) 
        sensor_map["gz_PMT9p"] = list(gz[500:525]) 
        sensor_map["gx_PMT10"] = list(gx[575:600]) 
        sensor_map["gy_PMT10"] = list(gy[575:600]) 
        sensor_map["gz_PMT10"] = list(gz[575:600]) 
        sensor_map["gx_PMT10p"] = list(gx[625:650]) 
        sensor_map["gy_PMT10p"] = list(gy[625:650]) 
        sensor_map["gz_PMT10p"] = list(gz[625:650]) 
        sensor_map["gx_PMT11"] = list(gx[700:725]) 
        sensor_map["gy_PMT11"] = list(gy[700:725]) 
        sensor_map["gz_PMT11"] = list(gz[700:725]) 
        sensor_map["gx_PMT11p"] = list(gx[750:775]) 
        sensor_map["gy_PMT11p"] = list(gy[750:775]) 
        sensor_map["gz_PMT11p"] = list(gz[750:775]) 
        sensor_map["gx_PMT4"] = list(gx[800:825]) 
        sensor_map["gy_PMT4"] = list(gy[800:825]) 
        sensor_map["gz_PMT4"] = list(gz[800:825]) 
        sensor_map["gx_PMT4p"] = list(gx[875:900]) 
        sensor_map["gy_PMT4p"] = list(gy[875:900]) 
        sensor_map["gz_PMT4p"] = list(gz[875:900]) 
        sensor_map["gx_PMT5"] = list(gx[925:950]) 
        sensor_map["gy_PMT5"] = list(gy[925:950]) 
        sensor_map["gz_PMT5"] = list(gz[925:950]) 
        sensor_map["gx_PMT5p"] = list(gx[975:1000]) 
        sensor_map["gy_PMT5p"] = list(gy[975:1000]) 
        sensor_map["gz_PMT5p"] = list(gz[975:1000]) 
        sensor_map["gx_PMT6"] = list(gx[1050:1075]) 
        sensor_map["gy_PMT6"] = list(gy[1050:1075]) 
        sensor_map["gz_PMT6"] = list(gz[1050:1075]) 
        sensor_map["gx_PMT6p"] = list(gx[1100:1125]) 
        sensor_map["gy_PMT6p"] = list(gy[1100:1125]) 
        sensor_map["gz_PMT6p"] = list(gz[1100:1125]) 
        sensor_map["gx_PMT7"] = list(gx[1150:1175]) 
        sensor_map["gy_PMT7"] = list(gy[1150:1175]) 
        sensor_map["gz_PMT7"] = list(gz[1150:1175]) 
        sensor_map["gx_PMT7p"] = list(gx[1225:1250]) 
        sensor_map["gy_PMT7p"] = list(gy[1225:1250]) 
        sensor_map["gz_PMT7p"] = list(gz[1225:1250]) 
        sensor_map["gx_PMT8"] = list(gx[1280:1305]) 
        sensor_map["gy_PMT8"] = list(gy[1280:1305]) 
        sensor_map["gz_PMT8"] = list(gz[1280:1305]) 
        # print(f"bmap {sensor_map}")
    elif run == 3:
        sensor_map["bx_PMT8"] = list(bx[0:25]) #PMT 8 Nort)h
        sensor_map["by_PMT8"] = list(by[200:225]) #PMT 8 Nort)h
        sensor_map["bz_PMT8"] = list(bz[200:225]) #PMT 8 Nort)h
        sensor_map["bx_PMT8p"] = list(bx[300:325]) #PMT 8 and 9 Nort)h
        sensor_map["by_PMT8p"] = list(by[300:325]) #PMT 8 and 9 Nort)h
        sensor_map["bz_PMT8p"] = list(bz[300:325]) 
        sensor_map["bx_PMT9"] = list(bx[400:425]) 
        sensor_map["by_PMT9"] = list(by[400:425]) 
        sensor_map["bz_PMT9"] = list(bz[400:425]) 
        sensor_map["bx_PMT9p"] = list(bx[500:525]) 
        sensor_map["by_PMT9p"] = list(by[500:525]) 
        sensor_map["bz_PMT9p"] = list(bz[500:525]) 
        sensor_map["bx_PMT10"] = list(bx[575:600]) 
        sensor_map["by_PMT10"] = list(by[575:600]) 
        sensor_map["bz_PMT10"] = list(bz[575:600]) 
        sensor_map["bx_PMT10p"] = list(bx[625:650]) 
        sensor_map["by_PMT10p"] = list(by[625:650]) 
        sensor_map["bz_PMT10p"] = list(bz[625:650]) 
        sensor_map["bx_PMT11"] = list(bx[700:725]) 
        sensor_map["by_PMT11"] = list(by[700:725]) 
        sensor_map["bz_PMT11"] = list(bz[700:725]) 
        sensor_map["bx_PMT11p"] = list(bx[750:775]) 
        sensor_map["by_PMT11p"] = list(by[750:775]) 
        sensor_map["bz_PMT11p"] = list(bz[750:775]) 
        sensor_map["bx_PMT4"] = list(bx[800:825]) 
        sensor_map["by_PMT4"] = list(by[800:825]) 
        sensor_map["bz_PMT4"] = list(bz[800:825]) 
        sensor_map["bx_PMT4p"] = list(bx[875:900]) 
        sensor_map["by_PMT4p"] = list(by[875:900]) 
        sensor_map["bz_PMT4p"] = list(bz[875:900]) 
        sensor_map["bx_PMT5"] = list(bx[925:950]) 
        sensor_map["by_PMT5"] = list(by[925:950]) 
        sensor_map["bz_PMT5"] = list(bz[925:950]) 
        sensor_map["bx_PMT5p"] = list(bx[975:1000]) 
        sensor_map["by_PMT5p"] = list(by[975:1000]) 
        sensor_map["bz_PMT5p"] = list(bz[975:1000]) 
        sensor_map["bx_PMT6"] = list(bx[1050:1075]) 
        sensor_map["by_PMT6"] = list(by[1050:1075]) 
        sensor_map["bz_PMT6"] = list(bz[1050:1075]) 
        sensor_map["bx_PMT6p"] = list(bx[1100:1125]) 
        sensor_map["by_PMT6p"] = list(by[1100:1125]) 
        sensor_map["bz_PMT6p"] = list(bz[1100:1125]) 
        sensor_map["bx_PMT7"] = list(bx[1150:1175]) 
        sensor_map["by_PMT7"] = list(by[1150:1175]) 
        sensor_map["bz_PMT7"] = list(bz[1150:1175]) 
        sensor_map["bx_PMT7p"] = list(bx[1225:1250]) 
        sensor_map["by_PMT7p"] = list(by[1225:1250]) 
        sensor_map["bz_PMT7p"] = list(bz[1225:1250]) 
        sensor_map["bx_PMT8"] = list(bx[1280:1305]) 
        sensor_map["by_PMT8"] = list(by[1280:1305]) 
        sensor_map["bz_PMT8"] = list(bz[1280:1305]) 
        # print(f"bmap {sensor_map}")
    file_name = f"./pole_rotation_measurement_run{run}.json"
    with open(file_name, "w") as file:
        json.dump(sensor_map, file, indent=4)

if __name__ == "__main__":
    main()