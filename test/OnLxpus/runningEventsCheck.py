#!/usr/bin/env python

import sys, os, string, re, copy
from optparse import OptionParser
from ROOT import *

# usage description
usage = """Usage: ./runningEventsCheck.py [options]\n
Example: ./runningEventsCheck.py -w LXBatch_Jobs \n
For more help: ./runningEventsCheck.py  --help
"""

# input parameters
parser = OptionParser(usage=usage)

parser.add_option("-w", "--main_workdir", dest="main_workdir", help="Main working directory", metavar="MAIN_WORKDIR")

(options, args) = parser.parse_args()

if not options.main_workdir:
    print usage
    sys.exit()

DataNames={"SingleElectron":"Evt_CutFlow_El", "SingleMu":"Evt_CutFlow_Mu"}
RunNames=["Run2012A-22Jan2013-v1_190645-193621",
          "Run2012B-22Jan2013-v1_193834-194833",
          "Run2012B-22Jan2013-v1_194834-195833",
          "Run2012B-22Jan2013-v1_195834-196531",
          "Run2012C-22Jan2013-v1_198049-199649",
          "Run2012C-22Jan2013-v1_199650-200649",
          "Run2012C-22Jan2013-v1_200650-201949",
          "Run2012C-22Jan2013-v1_201950-203002",
          "Run2012D-22Jan2013-v1_203709-205509",
          "Run2012D-22Jan2013-v1_205510-206509",
          "Run2012D-22Jan2013-v1_206510-207409",
          "Run2012D-22Jan2013-v1_207410-208686"
        ]
#HistName=["Evt_CutFlow_El","Evt_CutFlow_Mu"]
TreeName="SemiLeptanic"

for DataName in DataNames:
    totalEvts=0
    HistName=DataNames[DataName]
    print ">> [INFO] "+DataName
    for RunName in RunNames:
        rootname=options.main_workdir+"/"+DataName+"_"+RunName+".root"
        f = TFile(rootname)
        h = f.Get(TreeName+"/"+HistName)
        evts = h.GetBinContent(1)
        totalEvts+=evts
        print ">>        "+RunName+" "+str(int(evts))
    print ">>        Total "+str(int(totalEvts))

       



