#!/usr/bin/env python
import os, re, sys, shutil, commands
import math

dataPath='/BprimeKitNtuples/Production_CMSSW5311/CMSSW_5_3_11_data_8TeV_22Jan2013ReReco_AOD_v4'
#ref = open('singleElectron.txt', 'r')
ref = open('singleMuon.txt', 'r')
fileName='results'
for iref in ref:
	rfile = iref.strip()	
	print '>> %s...'%(rfile)
	script = open('./'+rfile+'.sh', 'w')
	listDir = commands.getoutput('ls '+dataPath+'/'+rfile+'/')
	files = listDir.split('\n')
	for file in files:
		if file.split('_')[0] == fileName:
			script.write('scp ntucms1.cern.ch:'+dataPath+'/'+rfile+'/'+file+' .\n')
			script.write('sleep 5\n')
	script.close()

