#!/usr/bin/env python

import sys, os, shutil, re, subprocess, time, glob
from optparse import OptionParser

rootName='SemiLeptanicAnalysis'
cmseos='/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select'

### * bsub Error description
errorMsg = { 
"Fail open root" : ["Fatal Root Error", "Failed to open the file"],
"EOS error"      : ["Error while doing the asyn writing", "error: target file was not created!"],
"Segmentation"   : ["Segmentation"],
"Over CPU/RAM"   : ["CPU time limit exceeded", "Aborted", "Killed", "bad_alloc"],
}

### * usage description
usage = """
>> [INFO] ./checkJobs.py [options]
>>        Example: ./checkJobs.py -w [workDir]
>>        Detial: . /checkJobs.py --help 
"""

### * Assistant functions
#
# 1. Check log and store the error according to 'errorMsg'
def checkError( outputDir, failJobDetail ):

    hasError=False
    for errType in errorMsg:
        failJobDetail[errType] = {} 
        for err in errorMsg[errType]:
            failJobDetail[errType][err] = [] 
            cmd = 'grep -r \''+err+'\' '+outputDir
            errInJobs = os.popen(cmd).read()
            if errInJobs != '':
                for errJob in errInJobs.split('\n'):
                    if errJob == '':
                        continue
                    job = re.search('\d+', re.search('job_\d+.',errJob).group(0)).group(0)
                    if int(job) not in failJobDetail[errType][err]:
                        failJobDetail[errType][err].append(int(job))
                        hasError=True
    return hasError


# 2. Check not done roots
def checkRoots( outputDir, nJobs, list_noDoneRoot, rootName ):

    originDir = os.getcwd()
    os.chdir(outputDir)
    roots = glob.glob( rootName+'_*.root' )

    i = 0
    while ( i < nJobs ):
        if rootName+'_'+str(i)+'.root' not in roots:
            list_noDoneRoot.append(i) 
        i+=1
    os.chdir(originDir)

    return len(list_noDoneRoot)


# 3. Check not done logs
def checkLogs( outputDir, nJobs, list_noDoneLogs ):

    originDir = os.getcwd()
    os.chdir(outputDir)
    logs = glob.glob( 'job_*.log' )
    i = 0
    while ( i < nJobs ):
        if 'job_'+str(i)+'.log' not in logs:
            list_noDoneLogs.append(i) 
        i+=1
    os.chdir(originDir)

    return len(list_noDoneLogs)


# 4. Fill infomations 
def storeInfo( dataPath, list_datasetDir, dict_failedInfo, list_noDoneRoot, list_noDoneLogs ):

    # Check if root output to eos
    isRootInEOS=False
    cmd = 'grep EOSPATH= '+dataPath+'/input/job_0.sh | awk -F \'"\' \'{print $2}\' | sed \'s/^\///g\''
    eospath=os.popen(cmd).read().strip()
    if eospath != '':
        print '>> ------------------------------------------------------- '
        ## * Mount eos if output in eos
        if not os.path.isdir('.eos') or os.listdir('.eos') == []:
            print '>> [INFO] Mounting eos to .eos ...' 
            cmd=cmseos+' -b fuse mount $PWD/.eos'
            os.system(cmd)
        else:
            print '>> [WARNING] .eos is not empty.' 
            print '>>           Did not mount eos again.'

        eospath='./.'+eospath
        isRootInEOS=True

    # Check root output name
    cmd = 'grep \'Copying file\' '+dataPath+'/input/job_0.sh | grep EOSPATH | awk \'{print $4}\''
    rootName = os.popen(cmd).read().strip().split('.root')[0]

    # Fill information
    list_datasetDir.append( len(glob.glob(dataPath+'/input/job_*.sh'))  )
    list_datasetDir.append( len(glob.glob(dataPath+'/output/job_*.log')))

    if isRootInEOS:
        list_datasetDir.append( len(glob.glob(eospath+'/'+rootName+'_*.root')) )
        checkRoots( eospath, list_datasetDir[0], list_noDoneRoot, rootName )
    else:
        list_datasetDir.append( len(glob.glob(dataPath+'/output/'+rootName+'_*.root')) )
        checkRoots( dataPath+'/output', list_datasetDir[0], list_noDoneRoot, rootName )

    checkLogs(  dataPath+'/output', list_datasetDir[0], list_noDoneLogs )
    checkError( dataPath+'/output', dict_failedInfo )

    # unMount eos if output in eos
    if isRootInEOS:
        if os.path.isdir('.eos') and os.listdir('.eos') != []:
            print '>> [INFO] unMounting eos from .eos ...' 
            cmd=cmseos+' -b fuse umount $PWD/.eos'
            os.system(cmd)
            print '>> ------------------------------------------------------- '
        print '>> [INFO] Output root in cmseos %s'%( '/'+eospath )
        print '>> ------------------------------------------------------- '

    return list_datasetDir[0]


# 5. Resubmit
def resubmit( path, job, queue ):

    doneResubmit=False
    if os.path.isdir(path) and 'output' in os.listdir(path) and 'input' in os.listdir(path):
        if os.path.isfile(path+'/output/job_'+str(job)+'.log'):
            cmd = 'mv '+path+'/output/job_'+str(job)+'.log '+path
            os.system(cmd)
        cmd = 'bsub -q '+queue + ' -o ' +path+'/output/job_'+str(job)+'.log'+' source '+path+'/input/job_'+str(job)+'.sh' 
        os.system(cmd)
        doneResubmit=True

    return doneResubmit


### * Main working func.
def main():

    ## * input parameters
    parser = OptionParser(usage=usage)
    parser.add_option("-w", "--workDir",     dest="workDir",     action='store',       help="Main working directory",       metavar="myWorkDir"                                   )
    parser.add_option("-d", "--dataset",     dest="dataset",     action='store',       help="Dataset in working directory", metavar="dataName",     default=''                    )
    parser.add_option("-q", "--queue",       dest="queue",       action='store',       help="LXBatch queue",                metavar="queueName",    default='cmscaf1nd'           )
    parser.add_option("-r", "--resubmit",    dest="resubmit",    action='store',       help="re-submit percific job",       metavar="resubmitJob",  default=None                  )
    parser.add_option("-R", "--resubmitAll", dest="resubmitAll", action='store_true',  help="re-submit all failed job",                             default=False                 )
    (options, args) = parser.parse_args()

    ## * make sure all necessary input parameters are provided
    if not (options.workDir):
        print usage
        sys.exit()

    allFiles=os.listdir(options.workDir)

    checkOneData=False
    doOneResubmit=False
    if options.resubmit and options.dataset=='':
        print '>> [ERROR] Please add dataset for resubmitting'
        print '>>        ./checkJob.py -w [workDir] -d [datasetName] -r [1 or 1,2,5...] (-q [queue])'
        sys.exit()
    elif options.resubmit:
        doOneResubmit=True

    if options.dataset != '':
        checkOneData=True
        if options.dataset in allFiles:
            f=options.workDir+'/'+options.dataset
            if not os.path.isdir(f) or 'output' not in os.listdir(f) or 'input' not in os.listdir(f):
                print '>> [ERROR] Not a work dataset dir : %s '%( options.dataset )
                print '>>        ./checkJob.py -w [workDir] -d [datasetName]'
                sys.exit()

    ## * Enter infomation
    resubmitJobs=[]
    if doOneResubmit :
        print '>> ------------------------------------------------------- '
        print '>> [INFO] Resubmitting...' 
        print '>>        Workspace : %s '%( options.workDir )
        print '>>        Data name : %s '%( options.dataset )
        print '>>        Job list  : %s '%( str(options.resubmit) )
        print '>> ------------------------------------------------------- '
        for i in options.resubmit.split(','):
            if i not in resubmitJobs:
                resubmitJobs.append(int(i))
    else:
        print '>> ------------------------------------------------------- '
        print '>> [INFO] Checking jot status...' 
        print '>>        Workspace : %s '%( options.workDir )
        if checkOneData:
            print '>>        Data name : %s '%( options.dataset )
        print '>> ------------------------------------------------------- '

    ## * Store the status of each datasets by log
    datasetDir={}
    noDoneRoot={}
    noDoneLogs={}
    failedInfo={}

    if doOneResubmit or checkOneData:
        # Only check one percific dataset
        fname=options.dataset
        datasetDir[fname] = [] 
        failedInfo[fname] = {}
        noDoneRoot[fname] = []
        noDoneLogs[fname] = []
        storeInfo( f, datasetDir[fname], failedInfo[fname], noDoneRoot[fname], noDoneLogs[fname] )
    else:
        # Check all datasets
        for fname in allFiles:
            f=options.workDir+'/'+fname
            if os.path.isdir(f) and 'output' in os.listdir(f) and 'input' in os.listdir(f):
                datasetDir[fname] = [] 
                failedInfo[fname] = {}
                noDoneRoot[fname] = []
                noDoneLogs[fname] = []
                storeInfo( f, datasetDir[fname], failedInfo[fname], noDoneRoot[fname], noDoneLogs[fname] )
    
    ## * Print and summerize all information
    nData=len(datasetDir)
    nDone=0

    # Go through each dataset
    for name in datasetDir:

        # 1. Count no done job ( including 'not found root', 'error in log' and 'not found log' )
        # 2. Simplize the info of error jobs
        noLogOrRoot=list(set(noDoneLogs[name]+noDoneRoot[name]))
        noYet=len(noLogOrRoot)
        sumErrJobs=[]
        failedJobs={}
        for errType in failedInfo[name]: # loop error type
            failedJobs[errType]=[]
            for err in failedInfo[name][errType]: # loop detail error msg
                # Simplize and store error jobs
                failedJobs[errType]=list(set(failedJobs[errType]+failedInfo[name][errType][err]))
                for fjob in failedInfo[name][errType][err]:
                    if fjob not in sumErrJobs:
                        sumErrJobs.append(int(fjob)) 
                        if fjob not in noLogOrRoot:  # Sum not done jobs
                            noYet+=1
                    
        # Print simple summary
        nJobs = int(datasetDir[name][0])
        nLogs = int(datasetDir[name][1])
        nRoot = int(datasetDir[name][2])
        nFail = len(sumErrJobs)
        nFinish = nJobs-noYet
        if nFinish == nJobs :
            print '> [DONE] %s '%( name )
            nDone+=1
        else:
            print '> %s '%( name )

        print '> N(Log)  : %d/%d '%( nLogs,   nJobs )
        print '> N(Root) : %d/%d '%( nRoot,   nJobs )
        print '> N(Done) : %d/%d '%( nFinish, nJobs )

        # Print error massage
        if nFail != 0:
            print '> N(Error): %d '%( nFail )
            print '  %-15s %s  '%( 'Sum fail jobs', str(sorted(sumErrJobs)).replace(" ", ""))
            for errType in failedJobs:
                if len(failedJobs[errType]) > 0:
                    print '  %-15s %s'%( errType,   str(sorted(failedJobs[errType])).replace(" ", ""))

        # Print not found logs 
        if len(noDoneLogs[name]) != 0:
            print '> %-15s %s'%( 'Not found logs',  str(sorted(noDoneLogs[name])).replace(" ", ""))

        # Print not found root 
        if len(noDoneRoot[name]) != 0:
            print '> %-15s %s'%( 'Not found roots', str(sorted(noDoneRoot[name])).replace(" ", ""))

        # Resubmit  
        if options.resubmitAll or doOneResubmit:
            if not doOneResubmit:
                resubmitJobs = sumErrJobs 
            print '> ReSubmiting all failed jobs %s '%( str(sorted(resubmitJobs)) )
            for job in resubmitJobs:
                if job in range(nJobs):
                    originDir = os.getcwd()
                    path = originDir+'/'+options.workDir+'/'+name
                    resubmit( path, job, options.queue )
                    print '  --- %d resubmitted!'%(job)         
                else:
                    print '  --- ERROR: %d out of [%d-%d]'%(job, 0, nJobs-1)         
                
        print '>> ------------------------------------------------------- '
    
    # Simple summary
    print '>> [INFO] Workspace  : %s   '%( options.workDir )
    print '>>        Root name  : %s   '%( rootName )
    if checkOneData:
        print '>>        Data name  : %s '%( options.dataset )
        print '>>        Done/Total : %d/1'%( nDone )
    else:
        print '>>        Done/Total : %d/%d'%( nDone, nData )

if __name__ == "__main__":
  main()

