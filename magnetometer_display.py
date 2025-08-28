#!/usr/bin/env python

#!/usr/bin/env python3

import time
import pickle as pkl
import csv

import tkinter as tk




import numpy as np

from iceboot.iceboot_session import getIcebootSession, startIcebootSession

# from optparse import OptionParser
# parser = OptionParser()
# parser.add_option("--host", dest="host", help="Ethernet host name or IP",
#                     default="localhost")
# parser.add_option("--port", dest="port", help="Ethernet port",
#                     default="5071")
# (options, args) = parser.parse_args()

# port = 5151 #mDOM m184 in Fieldhub01c4w2 WP0 5144+7
# port = 5146 #DEgg MB1 in Fieldhub01c4w2 WP0 5144+2
# port = 5148 #DEgg MB2 in Fieldhub01c4w2 WP0 5144+4

# port = 5156 #dEgg 2 in Fieldhub01c4w3 WP1 5152+4
# port = 5159 #mDOM DVT02 in Fieldhub01c4w3 WP1 5152+7

# port = 5066 # mDOM m091 in Fieldhub01c2w0 WP0 5064+2


# port = 5075 # DEgg 1 in Fieldhub01c2w0 WP0 5072+3


# host = options.host
# port = options.port

# device_map = {"5066":"m091","5151":"m184","5146":"DEggMB1","5148":"DEggMB2","5156":"DEgg2","5159":"mDOM_DVT_02","5075":"DEgg1"}
fl_chain_map = {"1":0x02,"2":0x04,"both":0x07}
led_map = {"1": 0x01,"2": 0x02,"3": 0x04,"4": 0x08,"5": 0x10,"all": 0x7F}

# timestr = time.strftime("%Y%m%d_%H%M%S")
timestr = time.strftime("%Y%m%d")
print(timestr)

def to_spherical(x,y,z):
    '''
    converts cartesian coordinates (x,y,z) to spherical (r,theta,phi)
    '''
    r = np.sqrt(x*x+y*y+z*z)
    theta = np.arctan2(np.sqrt(x*x + y*y),z)
    phi = np.arctan2(y,x)
    return r,theta,phi

root = tk.Tk()
display_label = tk.Label(root, text="Initial Value")
display_label.pack() # Or use .grid()
new_value = "Updated Value"
display_label.config(text=new_value)

def start_run():
    # session = getIcebootSession(host=host,port=port)
    session = getIcebootSession(devFile="/dev/tty.usbserial-FT5NADCZ", baudrate=1000000)
    # with open(f"../data/sensor_data/sensorsmDOMMB_West_{timestr}.txt","w") as text_file:
    with open(f"../data/sensor_data/sensorsmMB_{timestr}.txt","w") as text_file:
        # text_file.writelines(f"{'gx'} {'gy'} {'gz'} {'bx'} {'by'} {'bz'}\n")
        for i in range(0,10):
            try:
                bx,by,bz = session.readMagnetometerXYZ()
                r,theta,phi = to_spherical(bx,by,bz)
                print(f"Bx {bx*10**6:.1f}muT By {by*10**6:.1f}muT Bz {bz*10**6:.1f}muT ")
                print(f"B {r*10**6:.1f}muT theta {np.rad2deg(theta):.1f} phi {np.rad2deg(phi):.1f}")
            except ValueError:
                print(f"magnetometer reading in device mDOM MB fails")
                bx,by,bz = [np.nan,np.nan,np.nan]
            time.sleep(1)
            # print(f"{i} magnetometer readings {bx} {by} {bz}")
            # print(session.readAccelerometerTemperature())
            # try:
            #     gx,gy,gz = session.readAccelerometerXYZ()
            # except ValueError:
            #     print(f"accelerometer reading in device mDOM MB fails")
            #     gx,gy,gz = [np.nan,np.nan,np.nan]
            # print(f"{i} accelerometer readings {gx} {gy} {gz}")
            # text_file.writelines(f"{gx} {gy} {gz} {bx} {by} {bz}\n")

start_run()