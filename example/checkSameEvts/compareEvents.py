#!/usr/bin/env python
import os, re, sys, shutil
import math, ROOT
import numpy

RunNo = 0
LumiNo = 1
EvtNo = 2

fileNames = ["SingleMu_MuChannel", "SingleMu_ElChannel", "SingleElectron_MuChannel", "SingleElectron_ElChannel" ]
#fileNames = ["SingleElectron_MuChannel", "SingleMu_MuChannel_vs_SingleElectron_MuChannel" ]
#fileNames = ["SingleMu_MuChannel_vs_SingleElectron_MuChannel", "SingleElectron_MuChannel" ]
fileInfos = {}

print ">> [INFO] Recording events infomations"
for name in fileNames:
    print ">>        %s..."% name
    lines0 = open(name+".txt")
    lines1 = filter( None, ( line.strip() for line in lines0 ))                             # Remove empty line
    lines2 = filter( lambda x: not x.startswith('#'),  ( line.strip() for line in lines1 )) # Remove # comments lines
    lines3 = filter( None, ( line.split() for line in lines2 ))                             # Split RunNo, LumiNo, EvtNo
    fileInfos[name] = lines3
    size=len(fileInfos[name])
    i=1
    totalsame=0
    for info1 in fileInfos[name]:
        if  i==1 or i%5000==0:
            print ">>        %d/%d..."%( i, size) 

        same = 0
        for info2 in fileInfos[name]:
            if info1 == info2:
                same += 1
        if same > 1:
            totalsame += 1
            print ">> [INFO] %d count events: %s %s %s"%( same, info1[RunNo], info1[LumiNo], info1[EvtNo])
        i+=1
    print ">> [INFO] %d events double count"%( totalsame)

for i in range(0, len(fileNames)):
    for j in range( i+1, len(fileNames)):

        sameEvts=0
        #name1 = ""
        #name2 = ""
        #if len(fileInfos[fileNames[i]]) < len(fileInfos[fileNames[j]]):
        #    name1 = fileNames[i]
        #    name2 = fileNames[j]
        #else:
        #    name1 = fileNames[j]
        #    name2 = fileNames[i]

        name1 = fileNames[i]
        name2 = fileNames[j]
        info1 = fileInfos[name1]
        info2 = fileInfos[name2]

        outfile = open(name1+"_vs_"+name2+".txt", "w")
        outfile.write("# %6s %7s %12s\n"%( "RunNo", "LumiNo", "EvtNo" ))

        print ">> [INFO] Checking %25s:%-7d %s:%d"%( name1, len(info1), name2, len(info2) )
        for evt1 in info1:
            if evt1 in info2:
                sameEvts += 1
                outfile.write("%8s %7s %12s\n"%( evt1[RunNo], evt1[LumiNo], evt1[EvtNo] ))

        outfile.write("# %s:%d %s:%d\n"%( name1, len(info1), name2, len(info2) ))
        if sameEvts == 0:
            print ">>        No same events!"
            outfile.write("# No same events.\n")
        else:
            print ">>        %d same events."% sameEvts
            outfile.write("# %d same events\n"% sameEvts )

