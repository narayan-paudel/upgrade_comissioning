#!/usr/bin/env python

import os
import glob

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from pathlib import Path
home = str(Path.home())

from utils import to_spherical, tilt_angle, get_means, get_sub_df, to_spherical_list, to_spherical_list_rotated,meas_from_df
from calibrate_magnetometer_2d import corrected_ellipse

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams.update({'font.size': 10})

plotFolder = home+"/research_ua/icecube/upgrade/upgrade_comissioning/plots/"
# print(degg_list)

import pandas as pd


df_mDOM_117_2 = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/sensorsmDOMMB_20250908XNorth_Lab117.txt",header=0,sep=" ")
df_mDOM_Utah = pd.read_csv("/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/utah/mDOM_North_UtahScan_PMT8_16North.txt",header=0,sep=" ")


angles = [i*10 for i in range(0,37)]
colorsCustom = ['#8dd3c7','#bebada','#fb8072','#80b1d3','#fdb462','#b3de69']


def get_mean_B_rolling(df,angle):
    '''
    returns field in muT and degrees
    '''
    bx,by,bz = df[["bx","by","bz"]].values[angle:angle+10].T
    # print(f"bx {bx} by {by} bz {bz}")
    angle_list = []
    r,theta,phi = to_spherical_list(bx,by,bz)
    return np.mean(bx)*10**6,np.std(bx)*10**6,np.mean(by)*10**6,np.std(by)*10**6,np.mean(bz)*10**6,np.std(bz)*10**6,\
        np.mean(r)*10**6,np.std(r)*10**6,np.rad2deg(np.mean(theta)),np.rad2deg(np.std(theta)),np.rad2deg(np.mean(phi)),np.rad2deg(np.std(phi))

def extract_xyz_from_df(df):
    angles_list = []
    r_list = []
    r_std_list = []
    theta_list = []
    theta_std_list = []
    phi_list = []
    phi_std_list = []
    bx_list = []
    bx_std_list = []
    by_list = []
    by_std_list = []
    bz_list = []
    bz_std_list = []
    for iangle in angles:
        bx,bx_std,by,by_std,bz,bz_std,r,r_std,theta,theta_std,phi,phi_std = get_mean_B_rolling(df,iangle)
        angles_list.append(iangle)
        bx_list.append(bx)
        bx_std_list.append(bx_std)
        by_list.append(by)
        by_std_list.append(by_std)
        bz_list.append(bz)
        bz_std_list.append(bz_std)
        r_list.append(r)
        r_std_list.append(r_std)
        theta_list.append(theta)
        theta_std_list.append(theta_std)
        phi_list.append(phi)
        phi_std_list.append(phi_std)
    bx_list_new = []
    by_list_new = []
    bz_list_new = []
    angles_list_new = []
    for ibx,iby,ibz,iangle in zip(bx_list,by_list,bz_list,angles_list):
        if not np.isnan(ibx) and not np.isnan(iby) and not np.isnan(ibz):
            bx_list_new.append(ibx)
            by_list_new.append(iby)
            bz_list_new.append(ibz)
            angles_list_new.append(iangle)
    bx_list = bx_list_new
    by_list = by_list_new
    bz_list = bz_list_new
    return bx_list,by_list,bz_list







bx_list,by_list,bz_list = extract_xyz_from_df(df_mDOM_117_2)

points = np.vstack((bx_list, by_list, bz_list)).T
print(points.shape)
print(points.mean(axis=0))
centroid = points.mean(axis=0)











# ------------------------------------------------------------
# Utility: Fit plane using SVD
# ------------------------------------------------------------
def fit_plane(points):
    centroid = points.mean(axis=0)
    X = points - centroid
    _, _, Vt = np.linalg.svd(X)
    normal = Vt[2]
    basis1 = Vt[0]
    basis2 = Vt[1]
    return centroid, normal, basis1, basis2


# ------------------------------------------------------------
# Utility: algebraic ellipse fit in 2D (Fitzgibbon)
# ------------------------------------------------------------
def fit_ellipse_2d(x, y):
    D = np.vstack([x*x, x*y, y*y, x, y, np.ones_like(x)]).T
    S = np.dot(D.T, D)
    C = np.zeros([6, 6])
    C[0, 2] = C[2, 0] = 2
    C[1, 1] = -1

    eigvals, eigvecs = np.linalg.eig(np.linalg.inv(S) @ C)
    a = eigvecs[:, np.argmax(eigvals)]
    return a


# ------------------------------------------------------------
# Convert ellipse coefficients to parametric parameters
# ------------------------------------------------------------
def ellipse_params_from_coeffs(a):
    A, B, C, D, E, F = a

    # center
    M = np.array([[2*A, B],
                  [B, 2*C]])
    b = np.array([-D, -E])
    center = np.linalg.solve(M, b)
    x0, y0 = center

    # shifted constants
    F0 = (F + D*x0 + E*y0 +
          A*x0**2 + B*x0*y0 + C*y0**2)

    # rotation
    theta = 0.5 * np.arctan2(B, A - C)

    # semi-axes
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    Ao = A*cos_t**2 + B*cos_t*sin_t + C*sin_t**2
    Co = A*sin_t**2 - B*cos_t*sin_t + C*cos_t**2
    a_len = np.sqrt(-F0 / Ao)
    b_len = np.sqrt(-F0 / Co)

    return x0, y0, a_len, b_len, theta


# ------------------------------------------------------------
# RANSAC ellipse fitting in 3D
# ------------------------------------------------------------
def ransac_ellipse_3d(points, num_iter=2000, dist_thresh=0.05, min_inliers=30):
    best_inliers = []
    N = len(points)

    for _ in range(num_iter):
        # sample 6 points (minimum for ellipse)
        idx = np.random.choice(N, 6, replace=False)
        sample = points[idx]

        # 1) Fit plane
        centroid, normal, b1, b2 = fit_plane(sample)

        # 2) Project all points into plane
        rel = points - centroid
        px = rel @ b1
        py = rel @ b2

        # 3) Fit ellipse in 2D
        try:
            a = fit_ellipse_2d(px, py)
        except:
            continue

        # 4) Distance of points to ellipse (approximate algebraic distance)
        A, B, Cc, Dd, Ee, Ff = a
        dist = A*px*px + B*px*py + Cc*py*py + Dd*px + Ee*py + Ff
        dist = np.abs(dist)

        # 5) Count inliers
        inliers = np.where(dist < dist_thresh)[0]

        if len(inliers) > len(best_inliers):
            best_inliers = inliers

        if len(best_inliers) > min_inliers:
            break

    print(f"RANSAC inliers: {len(best_inliers)} / {N}")
    return best_inliers


# ------------------------------------------------------------
# Final fitting using inliers only
# ------------------------------------------------------------
def fit_ellipse_3d(points):
    # Best-fit plane
    centroid, normal, b1, b2 = fit_plane(points)

    # Project points into plane 2D
    rel = points - centroid
    px = rel @ b1
    py = rel @ b2

    # Fit 2D ellipse
    a = fit_ellipse_2d(px, py)
    x0, y0, A0, B0, theta = ellipse_params_from_coeffs(a)

    # Reconstruct ellipse in 3D
    t = np.linspace(0, 2*np.pi, 500)
    x2d = x0 + A0*np.cos(t)*np.cos(theta) - B0*np.sin(t)*np.sin(theta)
    y2d = y0 + A0*np.cos(t)*np.sin(theta) + B0*np.sin(t)*np.cos(theta)

    ellipse3d = centroid + np.outer(x2d, b1) + np.outer(y2d, b2)

    return ellipse3d, centroid, normal, (A0, B0), theta


ellipse3d, centroid, normal, (A0, B0), theta = fit_ellipse_3d(points)
print(f"Ellipse center: {centroid}")
print(f"Ellipse axes: {A0}, {B0}")
print(f"Ellipse normal: {normal}")  









# ============================================================
# Example (noisy + outliers)
# ============================================================
# true ellipse in 3D
# t = np.linspace(0, 2*np.pi, 500)
# a_true, b_true = 3.0, 1.2
# R = np.array([
#     [0.6, -0.3, 0.74],
#     [0.7,  0.7, 0.14],
#     [-0.4, 0.63, 0.66]
# ])
# ellipse_local = np.vstack((a_true*np.cos(t), b_true*np.sin(t), np.zeros_like(t)))
# pts_clean = (R @ ellipse_local).T
# pts_clean = points

# # Add noise
# pts_noisy = pts_clean + 0.05*np.random.randn(*pts_clean.shape)

# # Add outliers
# outliers = np.random.uniform(-4, 4, size=(40, 3))
# points = np.vstack([pts_noisy, outliers])

# Run RANSAC
# inliers = ransac_ellipse_3d(points)

# # Final refinement on inliers
# ellipse3d, ctr, normal, axes, theta = fit_ellipse_3d(points[inliers])


# # ============================================================
# # PLOT RESULT
# # ============================================================
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# ax.scatter(points[:,0], points[:,1], points[:,2], s=5, alpha=0.3)
# ax.scatter(points[inliers,0], points[inliers,1], points[inliers,2], s=10, color='r')
# ax.plot(ellipse3d[:,0], ellipse3d[:,1], ellipse3d[:,2], linewidth=3, color='k')

# ax.set_title("Robust RANSAC Ellipse Fit in 3D")
# plt.show()

def rotation_matrix_from_vectors(a, b):
    """
    Rotate vector a → vector b using Rodrigues' formula.
    """
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = np.dot(a, b)

    if np.isclose(c, 1.0):
        return np.eye(3)          # no rotation needed
    if np.isclose(c, -1.0):
        # opposite direction → 180° rotation
        # pick any perpendicular vector
        perp = np.array([1,0,0])
        if abs(a[0]) > 0.9:
            perp = np.array([0,1,0])
        v = np.cross(a, perp)
        v /= np.linalg.norm(v)
        K = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
        return np.eye(3) + 2*K@K   # 180° rotation

    # Rodrigues rotation formula
    s = np.linalg.norm(v)
    K = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])
    R = np.eye(3) + K + K @ K * ((1 - c) / (s**2))
    return R


# Rotate ellipse into horizontal plane
# ---------------------------------------------------------
def rotate_ellipse_to_horizontal(ellipse_pts,normal,centroid):
    # normal, centroid = fit_plane_normal(ellipse_pts)

    # We want the plane normal to become (0,0,1)
    target = np.array([0, 0, 1])
    R = rotation_matrix_from_vectors(normal, target)

    # Rotate around centroid
    pts_centered = ellipse_pts - centroid
    rotated = (R @ pts_centered.T).T + centroid

    return rotated, R, normal

rotated_ellipse, R_matrix, original_normal = rotate_ellipse_to_horizontal(points,normal,centroid)
ellipse3d_horizontal, centroid_horizontal, normal_horizontal, (A0, B0), theta_horizontal = fit_ellipse_3d(rotated_ellipse)


def plot_xyz_360(points,ellipse3d,centroid,normal,rotated_ellipse,ellipse3d_horizontal,centroid_horizontal,normal_horizontal, df_label, MB):
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0],projection='3d')
    ncolor = 0    
    ax.plot3D(points[:,0],points[:,1],points[:,2],"o",c=colorsCustom[ncolor],label=f"{df_label}",alpha=1)
    ax.plot(ellipse3d[:,0], ellipse3d[:,1], ellipse3d[:,2], linewidth=3, color=colorsCustom[ncolor])
    ax.plot3D(rotated_ellipse[:,0],rotated_ellipse[:,1],rotated_ellipse[:,2],"o",c=colorsCustom[ncolor+2],label=f"{df_label} rotated",alpha=1)
    ax.plot(ellipse3d_horizontal[:,0], ellipse3d_horizontal[:,1], ellipse3d_horizontal[:,2], linewidth=3, color=colorsCustom[ncolor+2])
    ax.quiver(centroid[0],centroid[1],centroid[2],normal[0], normal[1], normal[2], length=5, color=colorsCustom[ncolor], arrow_length_ratio=0.3)
    ax.quiver(centroid_horizontal[0],centroid_horizontal[1],centroid_horizontal[2],normal_horizontal[0], normal_horizontal[1], normal_horizontal[2], length=5, color=colorsCustom[ncolor+2], arrow_length_ratio=0.3)
    ax.plot3D([centroid[0]],[centroid[1]],[centroid[2]],"x",c="k",ms=10,mew=3,label="centroid")
    # Fit ellipse in 3D using RANSAC
    ax.tick_params(axis='both',which='both', direction='in', labelsize=16)
    ax.grid(True,alpha=0.6)
    ax.legend(loc="lower left",ncols=1,fontsize=16)
    ax.set_aspect('equal')
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" $B_x$ [$\mu$T]", fontsize=12)
    ax.set_ylabel(r" $B_y$ [$\mu$T]", fontsize=12)
    ax.set_zlabel(r" $B_z$ [$\mu$T]", fontsize=12)
    plt.savefig(plotFolder+f"/Bxyz_3d_ellipse{MB}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/Bxyz_3d_ellipse_{MB}.pdf",transparent=False,bbox_inches='tight')
    # plt.close()
    plt.show()



bx_list,by_list,bz_list = extract_xyz_from_df(df_mDOM_117_2)
# plot_xyz_360(bx_list,by_list,bz_list,"mDOM_117_2",MB="mDOM_117")
plot_xyz_360(points,rotated_ellipse,centroid,normal,rotated_ellipse,ellipse3d_horizontal,centroid_horizontal,normal_horizontal, "mDOM_117_2", MB="mDOM_117")

print(theta_horizontal,theta)
# print(np.anormal_horizontal,normal)
print(np.rad2deg(np.atan2(np.sqrt(normal_horizontal[1]**2 + normal_horizontal[0]**2),normal_horizontal[2])), np.rad2deg(np.atan2(np.sqrt(normal[1]**2 + normal[0]**2),normal[2])))


def extract_bx_by_bz(df):
    bx = df["bx"].values
    by = df["by"].values
    bz = df["bz"].values

    b_map = {}

    b_map["bx_step1"] = bx[:50] #PMT 8 North
    b_map["by_step1"] = by[:50] #PMT 8 North
    b_map["bz_step1"] = bz[:50] #PMT 8 North 
    b_map["bx_step2"] = bx[60:110] #PMT 9 North
    b_map["by_step2"] = by[60:110] #PMT 9 North
    b_map["bz_step2"] = bz[60:110] #PMT 9 North
    b_map["bx_step3"] = bx[150:200] #PMT 10 North
    b_map["by_step3"] = by[150:200] #PMT 10 North
    b_map["bz_step3"] = bz[150:200] #PMT 10 North
    b_map["bx_step4"] = bx[290:340] #PMT 11 North
    b_map["by_step4"] = by[290:340] #PMT 11 North
    b_map["bz_step4"] = bz[290:340] #PMT 11 North
    b_map["bx_step5"] = bx[380:430] #PMT 4 North
    b_map["by_step5"] = by[380:430] #PMT 4 North
    b_map["bz_step5"] = bz[380:430] #PMT 4 North
    b_map["bx_step6"] = bx[500:550] #PMT 5 North
    b_map["by_step6"] = by[500:550] #PMT 5 North
    b_map["bz_step6"] = bz[500:550] #PMT 5 North
    b_map["bx_step7"] = bx[650:700] #PMT 6 North
    b_map["by_step7"] = by[650:700] #PMT 6 North
    b_map["bz_step7"] = bz[650:700] #PMT 6 North
    b_map["bx_step8"] = bx[750:800] #PMT 7 North
    b_map["by_step8"] = by[750:800] #PMT 7 North
    b_map["bz_step8"] = bz[750:800] #PMT 7 North
    b_map["bx_step9"] = bx[860:910] #PMT 8 North
    b_map["by_step9"] = by[860:910] #PMT 8 North
    b_map["bz_step9"] = bz[860:910] #PMT 8 North
    b_map["bx_step10"] = bx[1000:1050] #PMT 9 North
    b_map["by_step10"] = by[1000:1050] #PMT 9 North
    b_map["bz_step10"] = bz[1000:1050] #PMT 9 North
    b_map["bx_step11"] = bx[1100:1150] #PMT 10 North
    b_map["by_step11"] = by[1100:1150] #PMT 10 North
    b_map["bz_step11"] = bz[1100:1150] #PMT 10 North
    b_map["bx_step12"] = bx[1200:1250] #PMT 11 North
    b_map["by_step12"] = by[1200:1250] #PMT 11 North
    b_map["bz_step12"] = bz[1200:1250] #PMT 11 North

    bx_mean_list = []
    by_mean_list = []
    bz_mean_list = []

    # step_list = [1,2,3,4,5,6,7,8,9,10,11,12,13]
    # step_list = [1,2,3,4,5,6,7,8,9,10,11,12]
    # step_list = [2,3,4,5,6,7,8,9,10]
    # step_list = [2,3,4,5,6,7,8,9,10] #first and last repeat to close circle
    step_list = [2,3,4,5,6,7,8,9] #first and last doesnot repeat to close circle Start from PMT 9
    # step_list = [1,2,3,4,5,6,7,8] #first and last doesnot repeat to close circle start from PMT 8

    for istep in step_list:
        bx_mean = np.mean(b_map[f"bx_step{istep}"])
        by_mean = np.mean(b_map[f"by_step{istep}"])
        bz_mean = np.mean(b_map[f"bz_step{istep}"])
        bx_mean_list.append(bx_mean*10**6)
        by_mean_list.append(by_mean*10**6)
        bz_mean_list.append(bz_mean*10**6)
    return bx_mean_list, by_mean_list, bz_mean_list


bx_mean_list_Utah, by_mean_list_Utah, bz_mean_list_Utah = extract_bx_by_bz(df_mDOM_Utah)
points_Utah = np.vstack((bx_mean_list_Utah, by_mean_list_Utah, bz_mean_list_Utah)).T
print(points_Utah.shape)

ellipse3d_utah, centroid_utah, normal_utah, (A0, B0), theta_utah = fit_ellipse_3d(points_Utah)

print(points_Utah)


rotated_ellipse_utah, R_matrix_utah, original_normal_utah = rotate_ellipse_to_horizontal(points_Utah,normal_utah,centroid_utah)
ellipse3d_horizontal_utah, centroid_horizontal_utah, normal_horizontal_utah, (A0, B0), theta_horizontal_utah = fit_ellipse_3d(rotated_ellipse_utah)

plot_xyz_360(points_Utah,ellipse3d_utah,centroid_utah,normal_utah,rotated_ellipse_utah,ellipse3d_horizontal_utah,centroid_horizontal_utah,normal_horizontal_utah,"mDOM_utah", MB="mDOM_Utah")
print(np.rad2deg(np.atan2(np.sqrt(normal_horizontal_utah[1]**2 + normal_horizontal_utah[0]**2),normal_horizontal_utah[2])), np.rad2deg(np.atan2(np.sqrt(normal_utah[1]**2 + normal_utah[0]**2),normal_utah[2])))

step_list = [2,3,4,5,6,7,8,9]
def plot_corrected_heading_360(bx_mean_list,by_mean_list,label=""):
    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(nrows=1,ncols=1, figure=fig)
    ax = fig.add_subplot(gs[0,0])
    # for df,dir in zip(df_list,dir_list):
    ncolor = 0
    bx_calibrated, by_calibrated = corrected_ellipse(bx_mean_list, by_mean_list)
    headings_original = [np.rad2deg(np.arctan2(by, bx)) for bx, by in zip(bx_mean_list, by_mean_list)]
    print(f"headings original: {headings_original}")
    headings_original = [i+360 if i<-1.0 else i for i in headings_original]
    print(f"headings original: {headings_original}")
    headings_corrected = [np.rad2deg(np.arctan2(byc, bxc)) for bxc,byc in zip(bx_calibrated, by_calibrated)]
    print(f"headings corrected: {headings_corrected}")
    headings_corrected = [i+360 if i<-1.0 else i for i in headings_corrected]
    # headings_corrected = [i+360 if 0<i<42 else i for i in headings_corrected]
    # headings_corrected = [i+360 if 0<i<34.0 else i for i in headings_corrected]
    # print(f"headings original: {headings_original}")
    print(f"headings corrected: {headings_corrected}")
    ###################################
    if label == "coarse":
        ax.plot(((np.asarray(step_list))-1)*45, headings_original, "-o", c=colorsCustom[ncolor], label=f"{'raw'}", alpha=1)
        ax.plot(((np.asarray(step_list))-1)*45, headings_corrected, "--o", c=colorsCustom[ncolor+2], label=f"{'calibrated'}", alpha=1)
        ax.set_xticks(np.linspace(0,360,9))
    else:
        ax.plot(range(len(headings_original)), headings_original, "-o", c=colorsCustom[ncolor], label=f"{'raw'}", alpha=1)
        ax.plot(range(len(headings_corrected)), headings_corrected, "--o", c=colorsCustom[ncolor+2], label=f"{'calibrated'}", alpha=1)
    #####################################
    ncolor += 1
    ax.tick_params(axis='both',which='both', direction='in', labelsize=20)
    ax.grid(True,alpha=0.6)
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xticks([str(int(i)) for i in np.linspace(0,360,9)])
    # ax.set_xlim(0,360)
    ###############################################
    ##############################################
    # ax.set_xticks(np.linspace(0,360,9)-roll)
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_xticks(np.linspace(0,360,9))
    # ax.set_yticks(np.linspace(0,360,9))
    # ax.set_xticks(np.linspace(0,540,13))
    # ax.set_yticks(np.linspace(-180,180,9))
    #############################################
    #############################################
    # ax.set_aspect('equal')
    ax.legend(ncols=2,fontsize=16)
    # ax.set_yticks(np.linspace(0,360,37))
    ax.set_xlabel(r" mDOM rotation [$^{\circ}$]", fontsize=20)
    ax.set_ylabel(r" $\phi$ [$^{\circ}$]", fontsize=20)
    plt.savefig(plotFolder+f"/orientation_with_Bxy_heading_Utah_{label}.png",transparent=False,bbox_inches='tight')
    plt.savefig(plotFolder+f"/orientation_with_Bxy_heading_Utah_{label}.pdf",transparent=False,bbox_inches='tight')
    plt.close()

print(f"ellipse data {rotated_ellipse_utah}")
plot_corrected_heading_360(rotated_ellipse_utah[:,0], rotated_ellipse_utah[:,1],label="coarse")