#!/usr/bin/env python
import os, re, sys, shutil, commands
import math

filePath='/BprimeKitNtuples/Production_CMSSW5311/CMSSW5_3_11_MC_AOD_v3'
ref = open('scp1.txt', 'r')
size=0
for rfile in ref:
	if rfile.strip() and not rfile.strip().startswith('#'):
		try:
			cmd1 = 'du '+filePath.strip()+'/'+rfile.strip()+'/'
			cmd2 = 'ls '+filePath.strip()+'/'+rfile.strip()+'/'

			d = commands.getoutput(cmd1).split('\n') 
			s = float(d[len(d)-1].split('\t')[0])
				
			ls = commands.getoutput(cmd2).split('\n')
			list=[]	
			for l in ls:
				if l.find('results') != -1 and l.find('root') != -1:
					list.append(l)
	
			#print str(len(list))+' files, '+str(s/1000)+'M '+rfile 
			#print '{0:<4d} files, {1:7.2f} M, {2:s}'.format(len(list), s/1000, rfile) 
			print '%-4d files, %7.2f M, %s'%( len(list), s/1000, rfile)
			size += s
		except:
			print '>> [ERROR] Not found %s'% rfile
sizeG = size/1000000
print str(sizeG)+'G in '+filePath

