#!/usr/bin/env python3

import numpy as np


def get_sub_df(df,angle):
    return df[df["angle"]==angle]

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
    #####to ensure phi is positive
    if phi < 0.0:
        phi_pos = phi + 2*np.pi
    else:
        phi_pos = phi
    return r,theta,phi_pos

def to_spherical_rotated(x,y,z):
    '''
    case when sensor is rotated with x axis pointing up (z direction)
    '''
    return to_spherical(z,y,x)

def to_spherical_list(x,y,z):
    '''
    converts list of cartesian coordinates (x,y,z) to spherical (r,theta,phi)
    '''
    r_list = []
    theta_list = []
    phi_list = []
    for ix,iy,iz in zip(x,y,z):
        r,theta,phi = to_spherical(ix,iy,iz)
        r_list.append(r)
        theta_list.append(theta)
        phi_list.append(phi)
    return r_list,theta_list,phi_list

def to_spherical_list_rotated(x,y,z):
    '''
    when sensor is rotated with x axis pointing up
    '''
    return to_spherical_list(z,y,x)

def get_mean(x,y,z):
    # print(imagnetometer_readings)
    ig_list = []
    itheta_list = []
    iphi_list = []
    for index, irow in iaccelerometer_readings.iterrows():
        # print(f"{index} {irow}")
        ig,itheta,iphi = to_spherical(irow["gx"],irow["gy"],irow["gz"])
        ig_list.append(ig)
        itheta_list.append(itheta)
        iphi_list.append(iphi)
        # print(f"gravity {ig} m/s2 {np.rad2deg(itheta)} {np.rad2deg(iphi)}")
    for ix in [ig_list,itheta_list,iphi_list]:
        accelerometer_sphe_dict[f"{iangle}"].append(np.mean(ix))
        accelerometer_sphe_dict[f"{iangle}"].append(np.std(ix))
    # print(f"{iaccelerometer_readings["gx"].values}")
    for jx in [iaccelerometer_readings["gx"].values,iaccelerometer_readings["gy"].values,iaccelerometer_readings["gz"].values]:
        accelerometer_cart_dict[f"{iangle}"].append(np.mean(jx))
        accelerometer_cart_dict[f"{iangle}"].append(np.std(jx))

def get_means(df):
    g_means = []
    g_stds = []
    gxy_means = []
    gxy_stds = []
    zenith_means = []
    zenith_stds = []
    azimuth_means = []
    azimuth_stds = []
    nx_means = []
    ny_means = []
    nz_means = []
    nx_stds = []
    ny_stds = []
    nz_stds = []
    # print(df.columns.values)
    for irot in df.columns.values[:]:
        # print(f"{irot} rotation")
        accelerometer_readings = df[[irot]].values    
        # print("accelerometer reading",accelerometer_readings)
        angles = []
        zenith = []
        azimuth = []
        g = []
        gxy = []
        nx_list = []
        ny_list = []
        nz_list = []
        for imeas in accelerometer_readings[:]:
            # print(imeas[0])
            # nx,ny,nz = imeas[0]
            print(imeas)
            nx,ny,nz = [float(jx.strip("[").strip("]")) for jx in imeas[0].split(",")]
            # print(tilt_angle(nx,ny,nz)*180/np.pi)
            angles.append(tilt_angle(nx,ny,nz)*180/np.pi)
            zenith.append(to_spherical(nx,ny,nz)[1]*180/np.pi)
            azimuth.append(to_spherical(nx,ny,nz)[2]*180/np.pi)
            g.append(to_spherical(nx,ny,nz)[0]*10**6)
            gxy.append(np.sqrt(nx**2+ny**2))
            nx_list.append(nx)
            ny_list.append(ny)
            nz_list.append(nz)
            # print(f"azimuth {azimuth}")
        # print(f"Rotation {irot} mean {np.mean(angles):.2f} std {np.std(angles):.2f}")s
        # print(f"Rotation {irot} mean zen {np.mean(zenith):.2f} std {np.std(zenith):.2f}"+
        #       f" mean azi {np.mean(azimuth):.2f} std {np.std(azimuth):.2f}"+
        #       f" mean B {np.mean(B):.2f} std {np.std(B):.2f}")
        azimuth = [i+360 if i<0 else i for i in azimuth]
        azimuth_means.append(np.mean(azimuth))
        azimuth_stds.append(np.std(azimuth))
        zenith_means.append(np.mean(zenith))
        zenith_stds.append(np.std(azimuth))
        g_means.append(np.mean(g))
        g_stds.append(np.std(g))
        gxy_means.append(np.mean(gxy))
        gxy_stds.append(np.std(gxy))
        nx_means.append(np.mean(nx_list))
        ny_means.append(np.mean(ny_list))
        nz_means.append(np.mean(nz_list))
        nx_stds.append(np.std(nx_list))
        ny_stds.append(np.std(ny_list))
        nz_stds.append(np.std(nz_list))

    # print(azimuth_means,g_means)
    return nx_means,nx_stds,ny_means,ny_stds,nz_means,nz_stds,g_means,g_stds,gxy_means,gxy_stds,azimuth_means,azimuth_stds, zenith_means,zenith_stds

def meas_from_df(df,measured_angles,sensor_rotation):
    # print(df)
    # print(df["gx"])
    magnetometer_cart_dict = {}
    magnetometer_sphe_dict = {}
    accelerometer_cart_dict = {}
    accelerometer_sphe_dict = {}

    for iangle in measured_angles:
        accelerometer_cart_dict[f"{iangle}"] = []
        accelerometer_sphe_dict[f"{iangle}"] = []
        magnetometer_cart_dict[f"{iangle}"] = []
        magnetometer_sphe_dict[f"{iangle}"] = []

    for iangle in measured_angles[:]:
        sub_df = get_sub_df(df,iangle)
        # print(iangle)
        # print(sub_df)
        bx,by,bz = sub_df[["bx","by","bz"]].values.T
        if not sensor_rotation:
            b,theta_b,phi_b = to_spherical_list(bx,by,bz)
        else:
            b,theta_b,phi_b = to_spherical_list_rotated(bx,by,bz)
        gx,gy,gz = sub_df[["gx","gy","gz"]].values.T
        if not sensor_rotation:
            g,theta_g,phi_g = to_spherical_list(gx,gy,gz)
        else:
            g,theta_g,phi_g = to_spherical_list_rotated(gx,gy,gz)
        for ix in [b,theta_b,phi_b]:
            magnetometer_sphe_dict[f"{iangle}"].append(np.mean(ix))
            magnetometer_sphe_dict[f"{iangle}"].append(np.std(ix))
        # print(f"{iaccelerometer_readings["gx"].values}")
        for jx in [bx,by,bz]:
            magnetometer_cart_dict[f"{iangle}"].append(np.mean(jx))
            magnetometer_cart_dict[f"{iangle}"].append(np.std(jx))
        for ix in [g,theta_g,phi_g]:
            accelerometer_sphe_dict[f"{iangle}"].append(np.mean(ix))
            accelerometer_sphe_dict[f"{iangle}"].append(np.std(ix))
        # print(f"{iaccelerometer_readings["gx"].values}")
        for jx in [gx,gy,gz]:
            accelerometer_cart_dict[f"{iangle}"].append(np.mean(jx))
            accelerometer_cart_dict[f"{iangle}"].append(np.std(jx))
    return accelerometer_cart_dict,accelerometer_sphe_dict,magnetometer_cart_dict,magnetometer_sphe_dict