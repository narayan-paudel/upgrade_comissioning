#!/usr/bin/env python3

import glob
import json
import re

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pathlib import Path
home = str(Path.home())

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']


plotFolder = home+"/research_ua/icecube/upgrade/upgrade_comissioning/plots/"
# print(degg_list)

geometry_folder = "/Users/epaudel/research_ua/icecube/software/upgrade_commissioning_scripts/geometry" #upgrade commissioning scripts repo
geometry_list = sorted(glob.glob(geometry_folder+"/string_??_geometry.json"))

upgrade_strings = [87,88,89,90,91,92]

print(geometry_list)

def get_geometry_xy(geometry_list):
    geometry_dict = {}
    for geometry_file in geometry_list[:]:
        with open(geometry_file, 'r') as f:
            geometry_data = json.load(f)
            string_num = int(re.search(r'string_(\d+)_geometry\.json', geometry_file).group(1))
            print(f"string number {string_num} geometry_data {geometry_data[0]['devices'][0]['x']}")
            x = geometry_data[0]['devices'][0]['x']
            y = geometry_data[0]['devices'][0]['y']

            geometry_dict[string_num] = [x,y]
    return geometry_dict

geometry_dict = get_geometry_xy(geometry_list)
print(f"geometry dict {geometry_dict}")

x_list = []
y_list = []
for istring in upgrade_strings[:]:
    x = geometry_dict[istring][0]
    y = geometry_dict[istring][1]
    x_list.append(x)
    y_list.append(y)


def get_distance(string1,string2):
    x1 = geometry_dict[string1][0]
    y1 = geometry_dict[string1][1]
    x2 = geometry_dict[string2][0]
    y2 = geometry_dict[string2][1]
    return np.sqrt((x2-x1)**2+(y2-y1)**2)

def get_angle(string1,string2,string3):
    x1 = geometry_dict[string1][0]
    y1 = geometry_dict[string1][1]
    x2 = geometry_dict[string2][0]
    y2 = geometry_dict[string2][1]
    x3 = geometry_dict[string3][0]
    y3 = geometry_dict[string3][1]
    v21 = np.array([x2-x1,y2-y1])
    v23 = np.array([x2-x3,y2-y3])
    cos_angle = np.dot(v21,v23)/(np.linalg.norm(v21)*np.linalg.norm(v23))
    angle_rad = np.arccos(cos_angle)
    angle_deg = np.degrees(angle_rad)
    return angle_deg

def plot_geometry(x_list,y_list,distances=None,angles=None):
    fig = plt.figure(figsize=(8,5))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    ax.scatter(x_list,y_list,label=f"geometry")

    ax.plot([x_list[0],x_list[1]],[y_list[0],y_list[1]],linestyle="--",color="gray",alpha=0.5)
    ax.plot([x_list[1],x_list[2]],[y_list[1],y_list[2]],linestyle="--",color="gray",alpha=0.5)
    ax.plot([x_list[2],x_list[3]],[y_list[2],y_list[3]],linestyle="--",color="gray",alpha=0.5)
    ax.plot([x_list[3],x_list[4]],[y_list[3],y_list[4]],linestyle="--",color="gray",alpha=0.5)
    ax.plot([x_list[4],x_list[5]],[y_list[4],y_list[5]],linestyle="--",color="gray",alpha=0.5)
    ax.plot([x_list[5],x_list[0]],[y_list[5],y_list[0]],linestyle="--",color="gray",alpha=0.5)
    ax.plot([x_list[0],x_list[2]],[y_list[0],y_list[2]],linestyle="--",color="gray",alpha=0.5)
    ax.plot([x_list[1],x_list[3]],[y_list[1],y_list[3]],linestyle="--",color="gray",alpha=0.5)
    ax.plot([x_list[1],x_list[4]],[y_list[1],y_list[4]],linestyle="--",color="gray",alpha=0.5)
    ax.plot([x_list[1],x_list[5]],[y_list[1],y_list[5]],linestyle="--",color="gray",alpha=0.5)
    if distances:
        dist_87_88 = get_distance(87,88)
        dist_88_89 = get_distance(88,89)
        dist_89_90 = get_distance(89,90)
        dist_90_91 = get_distance(90,91)
        dist_91_92 = get_distance(91,92)
        dist_92_87 = get_distance(92,87)
        dist_87_89 = get_distance(87,89)
        dist_88_90 = get_distance(88,90)
        dist_88_91 = get_distance(88,91)
        dist_88_92 = get_distance(88,92)
        ax.annotate(f"{dist_87_88:.1f}m", ((x_list[0]+x_list[1])/2, (y_list[0]+y_list[1])/2), textcoords="offset points", xytext=(0,0),color="green",fontsize="10", ha='center')
        ax.annotate(f"{dist_88_89:.1f}m", ((x_list[1]+x_list[2])/2, (y_list[1]+y_list[2])/2), textcoords="offset points", xytext=(0,0),color="green",fontsize="10", ha='center')
        ax.annotate(f"{dist_89_90:.1f}m", ((x_list[2]+x_list[3])/2, (y_list[2]+y_list[3])/2), textcoords="offset points", xytext=(0,0),color="green",fontsize="10", ha='center')   
        ax.annotate(f"{dist_90_91:.1f}m", ((x_list[3]+x_list[4])/2, (y_list[3]+y_list[4])/2), textcoords="offset points", xytext=(0,0),color="green",fontsize="10", ha='center')
        ax.annotate(f"{dist_91_92:.1f}m", ((x_list[4]+x_list[5])/2, (y_list[4]+y_list[5])/2), textcoords="offset points", xytext=(0,0),color="green",fontsize="10", ha='center')
        ax.annotate(f"{dist_92_87:.1f}m", ((x_list[5]+x_list[0])/2, (y_list[5]+y_list[0])/2), textcoords="offset points", xytext=(0,0),color="green",fontsize="10", ha='center')
        ax.annotate(f"{dist_87_89:.1f}m", ((x_list[0]+x_list[2])/2, (y_list[0]+y_list[2])/2), textcoords="offset points", xytext=(0,0),color="green",fontsize="10", ha='center')
        ax.annotate(f"{dist_88_90:.1f}m", ((x_list[1]+x_list[3])/2, (y_list[1]+y_list[3])/2), textcoords="offset points", xytext=(0,0),color="green",fontsize="10", ha='center')
        ax.annotate(f"{dist_88_91:.1f}m", ((x_list[1]+x_list[4])/2, (y_list[1]+y_list[4])/2), textcoords="offset points", xytext=(0,0),color="green",fontsize="10", ha='center')
        ax.annotate(f"{dist_88_92:.1f}m", ((x_list[1]+x_list[5])/2, (y_list[1]+y_list[5])/2), textcoords="offset points", xytext=(0,0),color="green",fontsize="10", ha='center')
    if angles:
        angle_87_88_89 = get_angle(87,88,89)
        angle_88_89_87 = get_angle(88,89,87)
        angle_89_87_88 = get_angle(89,87,88)
        angle_88_89_90 = get_angle(88,89,90)
        angle_89_90_88 = get_angle(89,90,88)
        angle_90_88_89 = get_angle(90,88,89)
        angle_88_90_91 = get_angle(88,90,91)
        angle_90_91_88 = get_angle(90,91,88)
        angle_91_88_90 = get_angle(91,88,90)
        angle_88_91_92 = get_angle(88,91,92)
        angle_91_92_88 = get_angle(91,92,88)
        angle_92_88_91 = get_angle(92,88,91)
        angle_88_92_87 = get_angle(88,92,87)
        angle_92_87_88 = get_angle(92,87,88)
        angle_87_88_92 = get_angle(87,88,92)

        ax.annotate(f"{angle_87_88_89:.0f}°", (x_list[1],y_list[1]), textcoords="offset points", xytext=(-20,-8), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_88_89_87:.0f}°", (x_list[2],y_list[2]), textcoords="offset points", xytext=(18,20), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_89_87_88:.0f}°", (x_list[0],y_list[0]), textcoords="offset points", xytext=(10,-15), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_88_89_90:.0f}°", (x_list[2],y_list[2]), textcoords="offset points", xytext=(30,8), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_89_90_88:.0f}°", (x_list[3],y_list[3]), textcoords="offset points", xytext=(-18,5), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_90_88_89:.0f}°", (x_list[1],y_list[1]), textcoords="offset points", xytext=(-5,-20), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_88_90_91:.0f}°", (x_list[3],y_list[3]), textcoords="offset points", xytext=(5,20), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_90_91_88:.0f}°", (x_list[4],y_list[4]), textcoords="offset points", xytext=(-20,-7), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_91_88_90:.0f}°", (x_list[1],y_list[1]), textcoords="offset points", xytext=(15,-9), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_88_91_92:.0f}°", (x_list[4],y_list[4]), textcoords="offset points", xytext=(-25,5), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_91_92_88:.0f}°", (x_list[5],y_list[5]), textcoords="offset points", xytext=(3,-20), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_92_88_91:.0f}°", (x_list[1],y_list[1]), textcoords="offset points", xytext=(18,2), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_88_92_87:.0f}°", (x_list[5],y_list[5]), textcoords="offset points", xytext=(-20,-15), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_92_87_88:.0f}°", (x_list[0],y_list[0]), textcoords="offset points", xytext=(30,0), ha='center',color="blue",fontsize=10)
        ax.annotate(f"{angle_87_88_92:.0f}°", (x_list[1],y_list[1]), textcoords="offset points", xytext=(-10,20), ha='center',color="blue",fontsize=10)

    for i,string in enumerate(upgrade_strings[:]):
        ax.annotate(f"{string}", (x_list[i],y_list[i]), textcoords="offset points", xytext=(0,10), ha='center')
    
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.4)
    # ax.legend(,ncols=3,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlim(10,85)
    ax.set_ylim(-95,-20)
    ax.set_aspect("equal")
    ax.set_xlabel(r"x [m]", fontsize=20)
    ax.set_ylabel(f"y [m]", fontsize=20)
    plt.savefig(plotFolder+f"/upgrade_geometry_dist_{distances}_angle_{angles}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/upgrade_geometry_dist_{distances}_angle_{angles}.pdf",transparent=False,bbox_inches='tight')
    plt.close()
plot_geometry(x_list,y_list,distances=True,angles=True)

# plot_geometry(x_list,y_list,distances=False,angles=False)
