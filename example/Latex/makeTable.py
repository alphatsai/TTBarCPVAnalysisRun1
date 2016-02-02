#!/usr/bin/env python
import os, re, sys, shutil
import math, ROOT
import numpy

from optparse import OptionParser

# usage description
usage = """
 Usage: """+sys.argv[0]+""" [options]
 Example: """+sys.argv[0]+""" --write results --input delphes.root --channel LepJets
 For more help: """+sys.argv[0]+""" --help 
"""

# Option parameters
parser = OptionParser(usage=usage)
parser.add_option("-i", "--input", dest="inputFile", 
                  help="Input file")
(options, args) = parser.parse_args()
if not options.inputFile:
    print usage
    sys.exit()
if not os.path.isfile(options.inputFile):
    print ">> [ERROR] Can't find "+options.inputFile+", or it's not a file..."
    sys.exit()

print ">> [INFO] Loaded %s"% options.inputFile
lines0 = open(options.inputFile)
lines1 = filter( None, ( line.strip() for line in lines0 ))                             # Remove empty line
lines2 = filter( lambda x: not x.startswith('#'),  ( line.strip() for line in lines1 )) # Remove # comments lines

nLines=len(lines2)
obsEl={"O2":"0", "O3":"0", "O4":"0", "O7":"0" }
systEl={"O2":{"+1Sigma":[],"-1Sigma":[]}, 
        "O3":{"+1Sigma":[],"-1Sigma":[]}, 
        "O4":{"+1Sigma":[],"-1Sigma":[]}, 
        "O7":{"+1Sigma":[],"-1Sigma":[]}}
sumSystEl={"O2":{"+1Sigma":"","-1Sigma":""}, 
        "O3":{"+1Sigma":"","-1Sigma":""}, 
        "O4":{"+1Sigma":"","-1Sigma":""}, 
        "O7":{"+1Sigma":"","-1Sigma":""}}
sumEl={"O2":{"+1Sigma":"","-1Sigma":""}, 
        "O3":{"+1Sigma":"","-1Sigma":""}, 
        "O4":{"+1Sigma":"","-1Sigma":""}, 
        "O7":{"+1Sigma":"","-1Sigma":""}}
obsMu={"O2":"0", "O3":"0", "O4":"0", "O7":"0" }
systMu={"O2":{"+1Sigma":[],"-1Sigma":[]}, 
        "O3":{"+1Sigma":[],"-1Sigma":[]}, 
        "O4":{"+1Sigma":[],"-1Sigma":[]}, 
        "O7":{"+1Sigma":[],"-1Sigma":[]}}
sumSystMu={"O2":{"+1Sigma":"","-1Sigma":""}, 
        "O3":{"+1Sigma":"","-1Sigma":""}, 
        "O4":{"+1Sigma":"","-1Sigma":""}, 
        "O7":{"+1Sigma":"","-1Sigma":""}}
sumMu={"O2":{"+1Sigma":"","-1Sigma":""}, 
        "O3":{"+1Sigma":"","-1Sigma":""}, 
        "O4":{"+1Sigma":"","-1Sigma":""}, 
        "O7":{"+1Sigma":"","-1Sigma":""}}
obsCo={"O2":"0", "O3":"0", "O4":"0", "O7":"0" }
systCo={"O2":{"+1Sigma":[],"-1Sigma":[]}, 
        "O3":{"+1Sigma":[],"-1Sigma":[]}, 
        "O4":{"+1Sigma":[],"-1Sigma":[]}, 
        "O7":{"+1Sigma":[],"-1Sigma":[]}}
sumSystCo={"O2":{"+1Sigma":"","-1Sigma":""}, 
        "O3":{"+1Sigma":"","-1Sigma":""}, 
        "O4":{"+1Sigma":"","-1Sigma":""}, 
        "O7":{"+1Sigma":"","-1Sigma":""}}
sumCo={"O2":{"+1Sigma":"","-1Sigma":""}, 
        "O3":{"+1Sigma":"","-1Sigma":""}, 
        "O4":{"+1Sigma":"","-1Sigma":""}, 
        "O7":{"+1Sigma":"","-1Sigma":""}}

latexObs={"O2":"\Otwo", "O3":"\Othree", "O4":"\Ofour", "O7":"\Oseven"}
latexSig={"+1Sigma":"+1", "-1Sigma":"-1"}

def outputElInfo( obsName, sigmaName ):
    if sigmaName.find("+1") > -1:
        print  "\multirow{2}{*}{%s} & \multirow{2}{*}{%s} & $%s \sigma$  & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s && %s & %s \\\\"%( 
        latexObs[obsName], obsEl[obsName], latexSig[sigmaName], systEl[obsName][sigmaName][1], systEl[obsName][sigmaName][2],systEl[obsName][sigmaName][3],systEl[obsName][sigmaName][4],systEl[obsName][sigmaName][5],systEl[obsName][sigmaName][6],systEl[obsName][sigmaName][7],systEl[obsName][sigmaName][8],systEl[obsName][sigmaName][9],systEl[obsName][sigmaName][10],systEl[obsName][sigmaName][11], sumSystEl[obsName][sigmaName], sumEl[obsName][sigmaName] 
        )
    else:
        print  " & & $%s \sigma$  & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s && %s & %s \\\\"%( 
        latexSig[sigmaName], systEl[obsName][sigmaName][1], systEl[obsName][sigmaName][2],systEl[obsName][sigmaName][3],systEl[obsName][sigmaName][4],systEl[obsName][sigmaName][5],systEl[obsName][sigmaName][6],systEl[obsName][sigmaName][7],systEl[obsName][sigmaName][8],systEl[obsName][sigmaName][9],systEl[obsName][sigmaName][10],systEl[obsName][sigmaName][11], sumSystEl[obsName][sigmaName], sumEl[obsName][sigmaName] 
        )
def outputMuInfo( obsName, sigmaName ):
    if sigmaName.find("+1") > -1:
        print  "\multirow{2}{*}{%s} & \multirow{2}{*}{%s} & $%s \sigma$  & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s && %s & %s \\\\"%( 
        latexObs[obsName], obsMu[obsName], latexSig[sigmaName], systMu[obsName][sigmaName][1], systMu[obsName][sigmaName][2],systMu[obsName][sigmaName][3],systMu[obsName][sigmaName][4],systMu[obsName][sigmaName][5],systMu[obsName][sigmaName][6],systMu[obsName][sigmaName][7],systMu[obsName][sigmaName][8],systMu[obsName][sigmaName][9],systMu[obsName][sigmaName][10],systMu[obsName][sigmaName][11], sumSystMu[obsName][sigmaName], sumMu[obsName][sigmaName] 
        )
    else:
        print  " & & $%s \sigma$  & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s && %s & %s \\\\"%( 
        latexSig[sigmaName], systMu[obsName][sigmaName][1], systMu[obsName][sigmaName][2],systMu[obsName][sigmaName][3],systMu[obsName][sigmaName][4],systMu[obsName][sigmaName][5],systMu[obsName][sigmaName][6],systMu[obsName][sigmaName][7],systMu[obsName][sigmaName][8],systMu[obsName][sigmaName][9],systMu[obsName][sigmaName][10],systMu[obsName][sigmaName][11], sumSystMu[obsName][sigmaName], sumMu[obsName][sigmaName] 
        )
def outputCoInfo( obsName, sigmaName ):
    if sigmaName.find("+1") > -1:
        print  "\multirow{2}{*}{%s} & \multirow{2}{*}{%s} & $%s \sigma$  & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s && %s & %s \\\\"%( 
        latexObs[obsName], obsCo[obsName], latexSig[sigmaName], systCo[obsName][sigmaName][1], systCo[obsName][sigmaName][2],systCo[obsName][sigmaName][3],systCo[obsName][sigmaName][4],systCo[obsName][sigmaName][5],systCo[obsName][sigmaName][6],systCo[obsName][sigmaName][7],systCo[obsName][sigmaName][8],systCo[obsName][sigmaName][9],systCo[obsName][sigmaName][10],systCo[obsName][sigmaName][11], sumSystCo[obsName][sigmaName], sumCo[obsName][sigmaName] 
        )
    else:
        print  " & & $%s \sigma$  & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s && %s & %s \\\\"%( 
        latexSig[sigmaName], systCo[obsName][sigmaName][1], systCo[obsName][sigmaName][2],systCo[obsName][sigmaName][3],systCo[obsName][sigmaName][4],systCo[obsName][sigmaName][5],systCo[obsName][sigmaName][6],systCo[obsName][sigmaName][7],systCo[obsName][sigmaName][8],systCo[obsName][sigmaName][9],systCo[obsName][sigmaName][10],systCo[obsName][sigmaName][11], sumSystCo[obsName][sigmaName], sumCo[obsName][sigmaName] 
        )

i=0
for line in lines2:
    if line.find('Electron channel') > -1 :
        obsName=lines2[i+1].split("(")[1].split(")")[0] #"Nominal Acp(O2) = -0.36"
        obsEl[obsName]=lines2[i+1].split()[3]
        systEl[obsName]["+1Sigma"]=lines2[i+3].split()
        systEl[obsName]["-1Sigma"]=lines2[i+4].split()
        sumSystEl[obsName]["+1Sigma"]=lines2[i+5].split()[2]
        sumSystEl[obsName]["-1Sigma"]=lines2[i+5].split()[3]
        sumEl[obsName]["+1Sigma"]=lines2[i+6].split()[2]
        sumEl[obsName]["-1Sigma"]=lines2[i+6].split()[3]
    elif line.find('Muon channel') > -1 :
        obsName=lines2[i+1].split("(")[1].split(")")[0] #"Nominal Acp(O2) = -0.36"
        obsMu[obsName]=lines2[i+1].split()[3]
        systMu[obsName]["+1Sigma"]=lines2[i+3].split()
        systMu[obsName]["-1Sigma"]=lines2[i+4].split()
        sumSystMu[obsName]["+1Sigma"]=lines2[i+5].split()[2]
        sumSystMu[obsName]["-1Sigma"]=lines2[i+5].split()[3]
        sumMu[obsName]["+1Sigma"]=lines2[i+6].split()[2]
        sumMu[obsName]["-1Sigma"]=lines2[i+6].split()[3]
    elif line.find('Combined channel') > -1 :
        obsName=lines2[i+1].split("(")[1].split(")")[0] #"Nominal Acp(O2) = -0.36"
        obsCo[obsName]=lines2[i+1].split()[3]
        systCo[obsName]["+1Sigma"]=lines2[i+3].split()
        systCo[obsName]["-1Sigma"]=lines2[i+4].split()
        sumSystCo[obsName]["+1Sigma"]=lines2[i+5].split()[2]
        sumSystCo[obsName]["-1Sigma"]=lines2[i+5].split()[3]
        sumCo[obsName]["+1Sigma"]=lines2[i+6].split()[2]
        sumCo[obsName]["-1Sigma"]=lines2[i+6].split()[3]
    i+=1   

#print systMu
 
print  "\multicolumn{17}{|c|}{Electron channel} \\\\"
print  "\hline"
outputElInfo( "O2", "+1Sigma")
outputElInfo( "O2", "-1Sigma")
print  "\hline"
outputElInfo( "O3", "+1Sigma")
outputElInfo( "O3", "-1Sigma")
print  "\hline"
outputElInfo( "O4", "+1Sigma")
outputElInfo( "O4", "-1Sigma")
print  "\hline"
outputElInfo( "O7", "+1Sigma")
outputElInfo( "O7", "-1Sigma")
print  "\hline"
print  "\multicolumn{17}{|c|}{Muon channel} \\\\"
print  "\hline"
outputMuInfo( "O2", "+1Sigma")
outputMuInfo( "O2", "-1Sigma")
print  "\hline"
outputMuInfo( "O3", "+1Sigma")
outputMuInfo( "O3", "-1Sigma")
print  "\hline"
outputMuInfo( "O4", "+1Sigma")
outputMuInfo( "O4", "-1Sigma")
print  "\hline"
outputMuInfo( "O7", "+1Sigma")
outputMuInfo( "O7", "-1Sigma")
print  "\hline"
print  "\multicolumn{17}{|c|}{Combined channel} \\\\"
print  "\hline"
outputCoInfo( "O2", "+1Sigma")
outputCoInfo( "O2", "-1Sigma")
print  "\hline"
outputCoInfo( "O3", "+1Sigma")
outputCoInfo( "O3", "-1Sigma")
print  "\hline"
outputCoInfo( "O4", "+1Sigma")
outputCoInfo( "O4", "-1Sigma")
print  "\hline"
outputCoInfo( "O7", "+1Sigma")
outputCoInfo( "O7", "-1Sigma")
print  "\hline"
