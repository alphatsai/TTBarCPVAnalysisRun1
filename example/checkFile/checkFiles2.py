#!/usr/bin/env python
import os, re, sys, shutil, commands
import math

filePath='/afs/cern.ch/user/j/jtsai/eos/cms/store/user/jtsai/TTBarCPV/bprimekits/EightTeV/TTJets_SemiLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A_ext-v1'
ref = open('TTJets_SemiLeptMGDecays.txt', 'r')
#filePath='/afs/cern.ch/user/j/jtsai/eos/cms/store/user/jtsai/TTBarCPV/bprimekits/EightTeV/TTJets_HadronicMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1'
#ref = open('TTJets_HadronicMGDecays.txt', 'r')
#filePath='/afs/cern.ch/user/j/jtsai/eos/cms/store/user/jtsai/TTBarCPV/bprimekits/EightTeV/TTJets_FullLeptMGDecays_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v2'
#ref = open('TTJets_FullLeptMGDecays.txt', 'r')

listDir = commands.getoutput('ls '+filePath)
files = listDir.split('\n')
fileName='results'
fileNum=[]
lostNum=[]
print '>> Loading files...'
for file in files:
	if file.split('_')[0] == fileName:
		fileNum.append(file.split('_')[1])
print '>> Loaded %d files'% len(fileNum)
print '>> Loading reference files...'

for rfile in ref:
	if rfile.split('_')[0] == fileName:
		idx = rfile.split('_')[1]
		if idx not in fileNum:
			lostNum.append(int(idx))
			print '>> [ERROR] Lost file: '+fileName+'_'+idx+'_'+rfile.split('_')[2]+'_'+rfile.split('_')[3].split('.')[0]+'.root';	
if len(lostNum) == 0:
	print '>> [INFO] All exist!'
else:	
	print '>> [INFO] Lost %d files'%len(lostNum)
