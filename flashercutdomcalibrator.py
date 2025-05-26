#!/usr/bin/env python

from I3Tray import *
from icecube import icetray, dataio, WaveCalibrator, dataclasses
from icecube.WaveCalibrator import DOMCalBaselineModule
from os.path import expandvars;
import sys
import os
from glob import glob
from string import atof,atoi

load("libdataclasses")
load("libphys-services")
load("libdataio")
load("libflasher-fill")
load("libI3Db")
load("libicepick")
#load("libquick-verification")
load("libpayload-parsing")
load("libdaq-decode")
#load("libpfauxiliary")
load("libflasher-timing")
#load("libFeatureExtractor")
load("libaligned-waveform")
load("libDOMcalibrator")
load("libflasher-verif")
load("libDomTools")
#load("libNFE")
load("libwavedeform")

tray = I3Tray()

#tools = expandvars("$I3_TOOLS")

#infile=sys.argv[1]


inice_rawdata                  = "InIceRawData"              # the name of the "InIce raw" DOMlaunch series frame-object
icetop_rawdata                 = "IceTopRawData"             # the name of the "IceTop raw" DOMlaunch series frame-object
inice_beacon                   = "InIceBeaconHits"           # the name of the "InIce beacon hits" DOMlaunch series frame-object
icetop_beacon                  = "IceTopBeaconHits"          # the name of the "IceTop beacon hits" DOMlaunch series frame-object
special_rawdata                = "SpecialRawData"            # the name of the "special DOM" DOMlaunch series frame-object
inice_clean_launches_name      = "InIceCleanedDomLaunches"   # the name of the "hit cleaned" InIce DOMlaunch series frame-object
icetop_clean_launches_name     = "IceTopCleanedDomLaunches"  # the name of the "hit cleaned" IceTop DOMlaunch series frame-object
inice_clean_launches_name_tmp  = "InIceCleanedDomLaunchesTmp"   # the name of the "hit cleaned" InIce DOMlaunch series frame-object
icetop_clean_launches_name_tmp = "IceTopCleanedDomLaunchesTmp"  # the name of the "hit cleaned" IceTop DOMlaunch series frame-object



inice_analysis_launches_name  = inice_clean_launches_name   # the name of the InIce DOMLaunch series to use in the analysis
icetop_analysis_launches_name = icetop_clean_launches_name  # the name of the IceTop DOMLaunches series to use in the analysis
hit_series                    = "InitialHitSeriesReco"      # the name of the "reco hit series" frame-object
pulse_series                  = "InitialPulseSeriesReco"    # the name of the "reco pulse series" frame-object

i3_header_name                = "I3DAQEventHeader"          # name of IceCube event header
i3_trigger_name               = "I3DAQTriggerHierarchy"     # name of IceCube trigger 
i3_buffer_name                = "I3DAQData"                 # name of raw IceCube buffer data

twr_header_name               = "TWRDAQEventHeader"         # name of TWR event header
twr_trigger_name              = "TWRDAQTriggerHierarchy"    # name of TWR trigger
twr_buffer_name               = "TWRDAQData"                # name of raw TWR buffer data
twr_raw_name                  = "TWRRawData"                # name of raw TWR frame object
twr_cleaned_name              = "TWRRawDataCleaned"         # name of cleaned raw TWR frame object
twr_reco_name                 = "TWRPulseSeriesReco"        # name of TWR reco pulse series, unshifted
twr_reco_shifted_name         = "TWRPulseSeriesRecoshifted" # name of TWR reco pulse series, shifted into IceCube time frame
twr_timeshift_name            = "TWRRawTimeCorrectionData"  # name of TWR time correction frame object
combined_reco_name            = "CombinedPulses"            # name of combined IceCube + TWR reco pulse series


infile=sys.argv[1]
#outfile=sys.argv[2]
# don't know if this works
outfile=sys.argv[2]

flasherfile=sys.argv[3]
waveformfile=sys.argv[4]
fstr=int(sys.argv[5])
fdom=int(sys.argv[6])

dbserver="dbs2.icecube.wisc.edu"
username="www"
workspace = expandvars("$I3_SRC")
#tray.AddService("I3XMLOMKey2MBIDFactory","omkey2mbid")(
#    ("infile",workspace+"/phys-services/resources/mainboard_ids.xml")
#    )

tray.AddService("I3DbOMKey2MBIDFactory","omkey2mbid")(
    ("host",dbserver),
    ("username",username),
    ("database","I3OmDb")
    )  


tray.AddService("I3FlasherFillServiceFactory","flashfill")(
    ("Hostname",dbserver),
    ("Username",username),
    ("DatabaseName","I3OmDb")
    )


tray.AddModule("I3Reader","readerfactory")(
    ("Filename", infile),
#    ("SkipUnregistered",True)
#    ("OmitGeometry",True),
#    ("OmitCalibration",True),
#    ("OmitStatus",True)
    )

tray.AddService("I3DbGeometryServiceFactory","geometry")(
    ("host",dbserver),
    ("username",username),
    ("database","I3OmDb")
    )

tray.AddService("I3DbCalibrationServiceFactory","calibration")(
    ("host",dbserver),
    ("username",username),
    ("database","I3OmDb")
    )

tray.AddService("I3DbDetectorStatusServiceFactory","status")(
    ("host",dbserver),
    ("username",username),
    ("database","I3OmDb")
    )


tray.AddModule( "I3MetaSynth", "challengeassumptions",
         GeometryService= "I3GeometryService",
         CalibrationService= "I3CalibrationService",
         DetectorStatusService= "I3DetectorStatusService"
         )

tray.AddModule("QConverter", "qify")
tray.AddService("I3PayloadParsingEventDecoderFactory","i3eventdecode")(
  #  ("Year",2008),
    ("headerid",i3_header_name),
    ("triggerid",i3_trigger_name),
    ("specialdataid",special_rawdata),
    ("specialdataoms",[OMKey(0,91),OMKey(0,92)]),
    ("flasherdataid","Flasher"),
    ("CPUDataID","BeaconHits")
    )

tray.AddModule("I3FrameBufferDecode","i3decode")(
          ("BufferID",i3_buffer_name)
          )



tray.AddModule("I3FlasherFillModule","flashfillmodule")

def myfilter(frame):
        if (frame.Has("flasher")):
            flasher = frame.Get("flasher")
            flash_om=[f.flashing_om.om for f in flasher]
	    flash_string=[f.flashing_om.string for f in flasher]
	    if ((flash_om[0]==fdom) and (flash_string[0]==fstr)):
		return 1
	    else:
		return 0
        else:
            return 0


tray.AddModule(myfilter, "Mynewfilter")

tray.AddModule("I3LCCleaning","LCClean_inice")(
      ("InIceInput","InIceRawData"),
        ("InIceOutput","InIceCleanData"))

#tray.AddSegment(WaveCalibrator.DOMSimulatorCalibrator, 'wavecal',
#               Launches="InIceCleanData",
#               )

#tray.AddModule("I3WaveformSplitter", "atwd_splitter",
#        Input="CalibratedWaveforms",
#        SplitATWDChannels=True
#        )


tray.AddModule("I3DOMcalibrator","merge")(
    ("InputRawDataName","InIceCleanData"),
    ("OutputATWDDataName", "InIceCalibratedATWD"),  
    ("CorrectPedestalDroopDualTau", False),           
    ("CalibrationMode", 0),  
    ("OutputFADCDataName", "InIceCalibratedFADC"),                    
    ("ATWDSaturationLevel",900))

outfile1 = "ATWD"+outfile
outfile2 = "FADC"+outfile

flasherfile1 = "ATWD"+flasherfile
flasherfile2 = "FADC"+flasherfile

waveformfile1 = "ATWD"+waveformfile
waveformfile2 = "FADC"+waveformfile

########tray.AddModule("Dump","dump")

tray.AddModule("I3AlignWaveform","aligned_waveform1")(
    #("ATWDReadoutName","UncalibratedATWD"),
    ("InIceLaunches","InIceCleanData"),
    #("CalibratedATWD","InIceCalibratedATWD")
    #("CalibratedFADC","InIceCalibratedFADC")
    ("OutputATWDDataName","InIceCalibratedATWD"),
    ("OutputFilename",outfile1),
    ("FlasherFilename",flasherfile1),
    ("WaveformFilename",waveformfile1),
    #("SumOnly",True),
    ("DumpWaveformFlag",False)
    )


tray.AddModule("I3AlignWaveform","aligned_waveform2")(
    ("InIceLaunches","InIceCleanData"),
    ("OutputATWDDataName","InIceCalibratedFADC"),
    ("OutputFilename",outfile2),
    ("FlasherFilename",flasherfile2),
    ("WaveformFilename",waveformfile2),
    ("SkipFirstNEvents",1),
#    ("SumOnly",True),
    ("DumpWaveformFlag",False)
    )



tray.AddModule("TrashCan","trash")

tray.Execute()
tray.Finish()
