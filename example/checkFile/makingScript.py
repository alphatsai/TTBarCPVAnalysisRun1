#!/usr/bin/env python
import os, re, sys, shutil, commands
import math

dataPath='/BprimeKitNtuples/Production_CMSSW5311/CMSSW_5_3_11_data_8TeV_22Jan2013ReReco_AOD_v4'
ref = open('data.txt', 'r')
eosPath='~/eos/cms/store/user/jtsai/TTBarCPV/bprimekits/EightTeV'
copy = open('copy.sh', 'w')
fileName='results'
for iref in ref:
	rfile = iref.strip()
	outfile = rfile+'.sh'
	copy.write('cp '+outfile+' '+eosPath+'/'+rfile+'/\n')		
	print '>> %s...'%(rfile)
	script = open(outfile, 'w')
	listDir = commands.getoutput('ls '+dataPath+'/'+rfile+'/')
	files = listDir.split('\n')
	for file in files:
		if file.split('_')[0] == fileName:
			script.write('scp ntucms1.cern.ch:'+dataPath+'/'+rfile+'/'+file+' .\n')
			script.write('sleep 5\n')
	script.close()
copy.close()		
