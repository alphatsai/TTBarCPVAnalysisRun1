#!/usr/bin/python
import os, re, sys, shutil
import math, time, array
from ROOT import *

f1 = TFile("MuonEfficiencies_Run2012ReReco_53X.root")
histNames=['DATA_over_MC_Loose_pt_abseta','DATA_over_MC_Tight_pt_abseta']
absEta_range=['<0.9','0.9-1.2','1.2-2.1','2.1-2.4'];
for histName in histNames:
    nominal=[]
    sigmaPos=[]
    sigmaNeg=[]
    for absEta in absEta_range:
        sf_nominal=[]
        sf_sigmaPos=[]
        sf_sigmaNeg=[]
        hName = histName+absEta
        print '>> [INFO] Get %s...'%hName
        h = f1.Get(hName)
        for pt in range(0,10):
            x = Double(0.) 
            y = Double(0.)
            h.GetPoint(pt,x,y)
            xMin=x-h.GetErrorXlow(pt)
            xMax=x+h.GetErrorXhigh(pt)
            yUp =h.GetErrorYhigh(pt)
            yLow=h.GetErrorYlow(pt)  
            print '          %3.0f-%-3.0f: %f +/- %f %f'%(xMin,xMax,y,yUp,yLow)
            sf_nominal.append(y)
            sf_sigmaPos.append(yUp)
            sf_sigmaNeg.append(yLow)
        nominal.append(sf_nominal)
        sigmaPos.append(sf_sigmaPos)
        sigmaNeg.append(sf_sigmaNeg)
    #print nominal
    #print sigmaPos
    #print sigmaNeg
            
