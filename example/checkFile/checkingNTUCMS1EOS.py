#!/usr/bin/env python
import os, re, sys, shutil, commands
import math

refPath='/BprimeKitNtuples/Production_CMSSW5311/CMSSW_5_3_11_data_8TeV_22Jan2013ReReco_AOD_v4'
eosPath='/eos/cms/store/group/dpg_ecal/alca_ecalcalib/ESAlignment/run1/Alpha/bprimekits/EightTeV/data'
files = open('samples.txt', 'r')
copy = open('copy.sh', 'w')
fileName='results'
refRoots=[]
eosRoots=[]
for iref in files:
	rfile = iref.strip()
	outfile = rfile+'.sh'
	copy.write('cp '+outfile+' ~'+eosPath+'/'+rfile+'/\n')		
	print '>> %s...'%(rfile)

	script = open(outfile, 'w')

	refroots = commands.getoutput('ls '+refPath+'/'+rfile+'/ | grep '+fileName)
	refRoots = refroots.split('\n')
	
	eosroots = commands.getoutput('/afs/cern.ch/project/eos/installation/cms/bin/eos.select ls '+eosPath+'/'+rfile+'/ | grep '+fileName)
	eosRoots = eosroots.split('\n')

	for root in refRoots:
		if root not in eosRoots:
			print '>> Not found %s'%root
			script.write('scp ntucms1.cern.ch:'+refPath+'/'+rfile+'/'+root+' .\n')
			script.write('sleep 5\n')
	script.close()
copy.close()		
