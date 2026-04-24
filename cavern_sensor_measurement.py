#!/usr/bin/env python3

import time
import numpy as np

from iceboot.iceboot_session import getIcebootSession


upgrade_string = 92
host = "fieldhub92"
port_list = [5111,5108,5106,5103,5100,5098,5127,5124,5122,5119,5116,5114,5143,5140,
             5138,5135,5132,5130,5159,5156,5154,5151,5148,5146,5175]

timestr = time.strftime("%Y%m%d")

def run_sensor_measurement(string,host,port):
    session = getIcebootSession(host=host,port=port)
    with open(f"/home/verical/epaudel/sensor_data/sensorsS{string}_{host}_p{port}_{timestr}.txt","w") as text_file:
        text_file.writelines(f"{'gx'} {'gy'} {'gz'} {'bx'} {'by'} {'bz'}\n")
        for i in range(0,100):
            # time.sleep(0.5)
            try:
                bx,by,bz = session.readMagnetometerXYZ() #tesla
            except ValueError:
                print(f"magnetometer reading in device mDOM MB fails")
                bx,by,bz = [np.nan,np.nan,np.nan]
            try:
                gx,gy,gz = session.readAccelerometerXYZ()
            except ValueError:
                print(f"accelerometer reading in device mDOM MB fails")
                gx,gy,gz = [np.nan,np.nan,np.nan]
            text_file.writelines(f"{gx} {gy} {gz} {bx} {by} {bz}\n")

for port in port_list:
    run_sensor_measurement(upgrade_string,host,port)