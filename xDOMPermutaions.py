#!/usr/bin/env python

import numpy as np
import pandas as pd

data = "/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/"
cmd_file = data + "strings_cmd/cmd_lite.xlsx"

df = pd.read_excel(cmd_file,header=None,usecols="B:H",skiprows=21,nrows=90)
# print(df)
n_neighbours = 2
xDOM_permutations_upgrade = []
for istring in np.arange(0,len(df.columns)):
    xDOM_list = df.iloc[:,istring]
    xDOM_permutations_string = []
    for i,imodule in enumerate(xDOM_list):
        # print(i)
        if i < 1:
            xDOM_above = "None"
        else:
            xDOM_above = xDOM_list[i-1]
        if i == len(xDOM_list)-1:
            xDOM_below = "None"
        else:
            xDOM_below = xDOM_list[i+1]
        # print(f"above {xDOM_above} this {imodule} below {xDOM_below}")
        if n_neighbours == 3:
            xDOM_permutations_string.append(f"{xDOM_above}_{imodule}_{xDOM_below}")
        elif n_neighbours == 2:
            xDOM_permutations_string.append(f"{imodule}_{xDOM_below}")
        else:
            print(f"n_neighbours is {n_neighbours} can only be 2 or 3")
    # print(xDOM_permutations_string)
    # print(f"For string\n {istring} {pd.Series(xDOM_permutations_string).value_counts()}")
    xDOM_permutations_upgrade += xDOM_permutations_string
print(f"For Upgrade\n {pd.Series(xDOM_permutations_upgrade).value_counts()}")
geoCal_permutations_upgrade = []
for icolumn in np.arange(0,len(df)):
    # print(f"this column print {icolumn}")
    this_column = df.iloc[icolumn].values
    geoCal_permutations_column = []
    for j,jxDOM in enumerate(this_column):
        if j < 1:
            xDOM_left = "None"
        else:
            xDOM_left = this_column[j-1]
        if j == len(this_column) -1:
            xDOM_right = "None"
        else:
            xDOM_right = this_column[j+1]
        if n_neighbours == 3:
            geoCal_permutations_column.append(f"{xDOM_left}_{jxDOM}_{xDOM_right}")
        elif n_neighbours == 2:
            geoCal_permutations_column.append(f"{jxDOM}_{xDOM_right}")
        else:
            print(f"n_neighbours is {n_neighbours} can only be 2 or 3")
        
    # print(f"For column \n {icolumn} {pd.Series(geoCal_permutations_column).value_counts()} ")
    geoCal_permutations_upgrade += geoCal_permutations_column
print(f"For Upgrade\n {pd.Series(geoCal_permutations_upgrade).value_counts()}")        
    
print("geometry calibration between strings")

