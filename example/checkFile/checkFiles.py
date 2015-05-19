#!/usr/bin/env python
import os, re, sys, shutil, commands
import math

filePath='/afs/cern.ch/user/j/jtsai/eos/cms/store/user/jtsai/TTBarCPV/bprimekits/EightTeV/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1'
listDir = commands.getoutput('ls '+filePath)
files = listDir.split('\n')
totalNum=512
fileName='results'
fileNum=[]
print '>> Loading files...'
for file in files:
	if file.split('_')[0] == fileName:
		fileNum.append(file.split('_')[1])
print '>> Checking files...'
for idx in range(1,512):
	if str(idx) not in fileNum:
		print '>> [ERROR] Lost file: '+str(idx);		
