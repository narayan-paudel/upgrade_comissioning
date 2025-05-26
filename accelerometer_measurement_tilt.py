#!/usr/bin/env python

import time
import numpy as np


from iceboot.iceboot_session import getIcebootSession
# session = getIcebootSession(devFile="/dev/tty.usbserial-FT5NADCZ", baudrate=1000000)




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
    return r,theta,phi

# gx,gy,gz = session.readAccelerometerXYZ()
# bx,by,bz = session.readMagnetometerXYZ()
# print(f"gx {gx} gy {gy} gz {gz}")
# titl_g = tilt_angle(gx,gy,gz)
# r_g,theta_g,phi_g = to_spherical(gx,gy,gz)
# print(f"bx {bx} by {by} bz {bz}")
# tilt_b = tilt_angle(bx,by,bz)
# r_b,theta_b,phi_b = to_spherical(bx,by,bz)
try:
    session = getIcebootSession(devFile="/dev/tty.usbserial-FT5NADCZ", baudrate=1000000)
except UnicodeDecodeError:
    print("UnicodeDecodeError: start session again")
    time.sleep(2)
else:
    print("no exception")
finally:                
    session = getIcebootSession(devFile="/dev/tty.usbserial-FT5NADCZ", baudrate=1000000)
# for irotation in range(0,361,10)[:]:
#     print(f"rotation {irotation}")
#     val = input("Enter y for next measurement: ")
#     if val == "y":
#         with open('/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/mDOMMB_accelerometer_magnetometerTiltMay16.txt', 'a+') as fh1:
#             for i in range(1,11):
#                 gx,gy,gz = session.readAccelerometerXYZ()
#                 bx,by,bz = session.readMagnetometerXYZ()
#                 fh1.write(f"{irotation}, {gx}, {gy}, {gz}, {bx}, {by}, {bz}\n")
#             input("lift tab A, input y")
#             for i in range(1,11):
#                 gx,gy,gz = session.readAccelerometerXYZ()
#                 bx,by,bz = session.readMagnetometerXYZ()
#                 fh1.write(f"{irotation}A, {gx}, {gy}, {gz}, {bx}, {by}, {bz}\n")
#             input("lift tab B, input y")
#             for i in range(1,11):
#                 gx,gy,gz = session.readAccelerometerXYZ()
#                 bx,by,bz = session.readMagnetometerXYZ()
#                 fh1.write(f"{irotation}B, {gx}, {gy}, {gz}, {bx}, {by}, {bz}\n")
#             input("lift tab C, input y")
#             for i in range(1,11):
#                 gx,gy,gz = session.readAccelerometerXYZ()
#                 bx,by,bz = session.readMagnetometerXYZ()
#                 fh1.write(f"{irotation}C, {gx}, {gy}, {gz}, {bx}, {by}, {bz}\n")
#             input("lift tab D, input y")
#             for i in range(1,11):
#                 gx,gy,gz = session.readAccelerometerXYZ()
#                 bx,by,bz = session.readMagnetometerXYZ()
#                 fh1.write(f"{irotation}D, {gx}, {gy}, {gz}, {bx}, {by}, {bz}\n")
# s = input("did you switch the main board, input y")
# if s == "y":
#     try:
#         session = getIcebootSession(devFile="/dev/tty.usbserial-FT5NADCZ", baudrate=1000000)
#     except UnicodeDecodeError:
#         print("UnicodeDecodeError: start session again")
#         time.sleep(2)
#     else:
#         print("no exception")
#     finally:                
#         session = getIcebootSession(devFile="/dev/tty.usbserial-FT5NADCZ", baudrate=1000000)
for irotation in range(0,361,10)[:]:
    print(f"rotation {irotation}")
    val = input("Enter y for next measurement: ")
    if val == "y":
        with open('/Users/epaudel/research_ua/icecube/upgrade/upgrade_comissioning/data/sensor_data/mMB_accelerometer_magnetometerTiltMay16.txt', 'a+') as fh2:
            for i in range(1,11):
                gx,gy,gz = session.readAccelerometerXYZ()
                bx,by,bz = session.readMagnetometerXYZ()
                fh2.write(f"{irotation} {gx} {gy} {gz} {bx} {by} {bz}\n")
            input("lift tab A, input y")
            for i in range(1,11):
                gx,gy,gz = session.readAccelerometerXYZ()
                bx,by,bz = session.readMagnetometerXYZ()
                fh2.write(f"{irotation}A {gx} {gy} {gz} {bx} {by} {bz}\n")
            input("lift tab B, input y")
            for i in range(1,11):
                gx,gy,gz = session.readAccelerometerXYZ()
                bx,by,bz = session.readMagnetometerXYZ()
                fh2.write(f"{irotation}B {gx} {gy} {gz} {bx} {by} {bz}\n")
            input("lift tab C, input y")
            for i in range(1,11):
                gx,gy,gz = session.readAccelerometerXYZ()
                bx,by,bz = session.readMagnetometerXYZ()
                fh2.write(f"{irotation}C {gx} {gy} {gz} {bx} {by} {bz}\n")
            input("lift tab D, input y")
            for i in range(1,11):
                gx,gy,gz = session.readAccelerometerXYZ()
                bx,by,bz = session.readMagnetometerXYZ()
                fh2.write(f"{irotation}D {gx} {gy} {gz} {bx} {by} {bz}\n")
   

