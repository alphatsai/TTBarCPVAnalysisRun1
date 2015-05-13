#!/usr/bin/env python

import sys, os, shutil, re, subprocess
from optparse import OptionParser


def make_filenamelist(input_dir):

    proc = subprocess.Popen( [ '/afs/cern.ch/project/eos/installation/cms/bin/eos.select', 'ls', input_dir ], stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
    output = proc.communicate()[0]
    if proc.returncode != 0:
        print output
        sys.exit(1)

    return output.splitlines()


def process_input_dir(input_dir, match, filelist):

    input_dir = input_dir.rstrip('/')+'/'
    filenamelist = make_filenamelist(input_dir)

    path = input_dir;
    name = ''
    jobdict = {}

    for filename in filenamelist:
        if( not re.search('.root$', filename) ):
            continue
        if ( match!=None and not re.search(match, filename) ):
            continue
        if re.search('_\d+_\d+_\w+.root', filename):
          m1 = re.search('_\d+_\d+_\w+.root', filename)
        elif re.search('_\d+.root', filename):
          m1 = re.search('_\d+.root', filename)
        if name=='':
          if re.search('_\d+_\d+_\w+.root', filename):
            name = re.split('_\d+_\d+_\w+.root', filename)[0]
          elif re.search('_\d+.root', filename):
            name = re.split('_\d+.root', filename)[0]
        jobstring = filename[m1.start():].lstrip('_').replace('.root','').split('_')
        job = int(jobstring[0])
        if len(jobstring)>1:
          if job not in jobdict.keys():
              jobdict[job] = []
              jobdict[job].append([int(jobstring[1])])
              jobdict[job].append([jobstring[2]])
          else:
              jobdict[job][0].append(int(jobstring[1]))
              jobdict[job][1].append(jobstring[2])
        else:
          if job not in jobdict.keys():
            jobdict[job] = []

    jobs = jobdict.keys()
    if( len(jobs)==0 ):
        print 'No matching .root files found'
        sys.exit()

    jobs.sort()
    for job in jobs:
        if len(jobdict[job])>1:
          maxsub = max(jobdict[job][0])
          filename = (path+name+'_%i_%i_%s.root')%(job, maxsub, jobdict[job][1][jobdict[job][0].index(maxsub)])
        else:
          filename = (path+name+'_%i.root')%(job)
        filelist.append(filename)

    return


cfi_template = """FileNames = [
INPUT_FILES
]
"""


bash_template = """#!/bin/bash

BATCHDIR=${PWD}
EOSPATH="EOS_PATH"

#export SCRAM_ARCH=slc5_amd64_gcc462
export SCRAM_ARCH=slc6_amd64_gcc481
cd MAIN_WORKDIR
eval `scram runtime -sh`

cp -v MAIN_WORKDIR/CMSSW_cfg.py $BATCHDIR/CMSSW_cfg.py
#cp -v MAIN_WORKDIR/inputMatrixElements_cfi.py $BATCHDIR/
cp -v DATASET_WORKDIR/input/inputFiles_JOB_NUMBER_cfi.py $BATCHDIR/inputFiles_cfi.py

cd $BATCHDIR
echo "Running CMSSW job"
cmsRun CMSSW_cfg.py CFG_PARAMETERS
exitcode=$?

if [ ! -z $EOSPATH -a $EOSPATH!=" " ]
then
  echo "Copying file  OUTPUT_FILENAME.root to " ${EOSPATH}
  source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh 
  /afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select cp OUTPUT_FILENAME.root ${EOSPATH}/OUTPUT_FILENAME_JOB_NUMBER.root
  #cmsStage -f OUTPUT_FILENAME.root 
else 
  echo "Copying file  OUTPUT_FILENAME.root to DATASET_WORKDIR/output/OUTPUT_FILENAME_JOB_NUMBER.root"
  cp -v OUTPUT_FILENAME.root DATASET_WORKDIR/output/OUTPUT_FILENAME_JOB_NUMBER.root
fi

exit $exitcode

"""


# usage description
usage = """Usage: ./createAndSubmitJobs.py [options]\n
Example: ./createAndSubmitJobs.py -w LXBatch_Jobs -d datasetList.txt -c btagvalidation_cfg.py\n
For more help: ./createAndSubmitJobs.py --help
"""

def main():
  # input parameters
  parser = OptionParser(usage=usage)

  parser.add_option("-w", "--main_workdir", dest="main_workdir", action='store', help="Main working directory", metavar="MAIN_WORKDIR")
  parser.add_option("-d", "--dataset_list", dest="dataset_list", action='store', help="Text file containing a list of datasets to be processed", metavar="DATASET_LIST")
  parser.add_option("-o", "--output_filename", dest="output_filename", action='store', default='ESSkim', help="Output ROOT filename (Default set to ESSkim)", metavar="OUTPUT_FILENAME")
  parser.add_option("-E", "--eos_path", dest="eos_path", action='store', help="EOS path to copy output files to (This parameter is optional)", metavar="EOS_PATH")
  parser.add_option('-m', '--match', dest="match", action='store', help='Only files containing the MATCH string in their names will be considered (This parameter is optional)', metavar='MATCH')
  parser.add_option("-c", "--cmssw_cfg", dest="cmssw_cfg", action='store', help="CMSSW configuration file", metavar="CMSSW_CFG")
  parser.add_option('-f', '--fraction', dest='fraction', action='store', default='1.0', help='Fraction of files to be processed. Default value is 1 (This parameter is optional)', metavar='FRACTION')
  parser.add_option("-q", "--queue", dest="queue", action='store', default='1nh', help="LXBatch queue (choose among cmst3 8nm 1nh 8nh 1nd 1nw). Default is '1nh' (This parameter is optional)", metavar="QUEUE")
  parser.add_option("-n", "--no_submission", dest="no_submission", action='store_true', default=False, help="Create the necessary configuration files and skip the job submission (This parameter is optional)")

  (options, args) = parser.parse_args()

  # make sure all necessary input parameters are provided
  if not (options.main_workdir and options.dataset_list and options.cmssw_cfg):
    print usage
    sys.exit()

  main_workdir = options.main_workdir
  dataset_list = options.dataset_list
  cmssw_cfg = options.cmssw_cfg

  # redefine main_workdir as an absolute path (if not defined in such form already)
  if not re.search("^/", main_workdir):
    main_workdir = os.path.join(os.getcwd(),main_workdir)

  # create the main working directory
  if not os.path.exists(main_workdir) : 
    os.mkdir(main_workdir)
  else :
    print 'Warning ' + main_workdir + ' exists. '

  # copy the dataset list file to the main_workdir
  shutil.copyfile(dataset_list,os.path.join(main_workdir,'datasetList.txt'))

  # copy the CMSSW cfg file to the cfg_files_dir
  shutil.copyfile(cmssw_cfg,os.path.join(main_workdir,'CMSSW_cfg.py'))

  # look for pileup distribution files and copy them into main_workdir
  cfg_dirname = os.path.dirname(cmssw_cfg)
  if cfg_dirname=='':
    cfg_dirname = os.getcwd()
  for filename in os.listdir(cfg_dirname):
    if not os.path.isfile(os.path.join(cfg_dirname,filename)):
      continue
    #if re.search("inputMatrixElements_cfi.py", filename):
    #  shutil.copy(os.path.join(cfg_dirname,filename),main_workdir)

  # open and read the dataset_list file
  dataset_list_file = open(dataset_list,"r")
  dataset_list_lines = dataset_list_file.readlines()

  # loop over datasets
  for line in dataset_list_lines:
    line_elements = line.split()
    if (len(line_elements)==0 or line_elements[0][0]=='#'): continue

    output_filename = options.output_filename
    cfg_parameters = ''
    if( len(line_elements)>3 ):
      cfg_parameters = line_elements[3]
      if( line_elements[3].split('=')[0]=='outFilename' ): output_filename = line_elements[3].split('=')[1].replace('.root','')
      for par in range(4,len(line_elements)):
        cfg_parameters = cfg_parameters + ' ' + line_elements[par]
        if( line_elements[par].split('=')[0]=='outFilename' ): output_filename = line_elements[par].split('=')[1].replace('.root','')

    dataset = line_elements[0].lstrip('/').replace('/','__')
    print 'Processing ' + line_elements[0]

    dataset_workdir = os.path.join(main_workdir,dataset)

    # create the dataset working directory
    os.mkdir(dataset_workdir)
    os.mkdir(os.path.join(dataset_workdir,'input'))
    os.mkdir(os.path.join(dataset_workdir,'output'))

    eos_path = ' '
    if ( options.eos_path ):
      eos_path = ''.join([options.eos_path,"/",dataset]) 
      proc = subprocess.Popen( [ '/afs/cern.ch/project/eos/installation/cms/bin/eos.select', 'mkdir', eos_path ], stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
      proc = subprocess.Popen( [ '/afs/cern.ch/project/eos/installation/cms/bin/eos.select', 'chmod 775', eos_path ], stdout = subprocess.PIPE, stderr = subprocess.STDOUT )

    filelist = []
    process_input_dir(line_elements[2], options.match, filelist)

    ##################################################
    njobs = line_elements[1]
    numfiles = int(len(filelist)*float(options.fraction))
    ijobmax=int(njobs)
    if ijobmax > numfiles:
       ijobmax=numfiles
       print '  Number of jobs requested exceeds the total number of input files.\n  The number of jobs set to ' + str(ijobmax) + ' to match the number of input files'
    filesperjob = int(numfiles/ijobmax)
    if numfiles%ijobmax!=0:
       filesperjob = filesperjob + 1
       ijobmax = int(numfiles/filesperjob)
       if numfiles%filesperjob!=0:
           ijobmax = ijobmax + 1
       if ijobmax != int(njobs):
           print '  Could not create the exact number of jobs requested.\n  For optimal job splitting, the number of jobs set to ' + str(ijobmax)
    #################################################
    ifile = 0
    for ijob in range(ijobmax):
       # prepare the list of input files
       input_files = '    \'root://eoscms.cern.ch/' + filelist[ifile] + '\''
       for i in range(filesperjob-1):
           if ifile>(numfiles-2):
               break
           ifile = ifile + 1
           input_files += (',\n    \'root://eoscms.cern.ch/' + filelist[ifile] + '\'')
       ifile = ifile + 1

       ## create input cfi file
       input_files_cfi = open(os.path.join(dataset_workdir,'input','inputFiles_' + str(ijob) + '_cfi.py'),'w')
       input_files_cfi.write(re.sub('INPUT_FILES',input_files,cfi_template))
       input_files_cfi.close()

       ## create Bash script
       bash_script = open(os.path.join(dataset_workdir,'input','job_' + str(ijob) + '.sh'),'w')
       bash_script_content = re.sub('MAIN_WORKDIR',main_workdir,bash_template)
       bash_script_content = re.sub('DATASET_WORKDIR',dataset_workdir,bash_script_content)
       bash_script_content = re.sub('JOB_NUMBER',str(ijob),bash_script_content)
       bash_script_content = re.sub('CFG_PARAMETERS',cfg_parameters,bash_script_content)
       bash_script_content = re.sub('OUTPUT_FILENAME',output_filename,bash_script_content)
       bash_script_content = re.sub('EOS_PATH',eos_path,bash_script_content)
       bash_script.write(bash_script_content)
       bash_script.close()

       if(not options.no_submission):
         cmd = 'bsub -q ' + options.queue + ' -o ' + os.path.join(dataset_workdir,'output','job_' + str(ijob) + '.log') + ' source ' + os.path.join(dataset_workdir,'input','job_' + str(ijob) + '.sh')
         os.system(cmd)

  # close all open files
  dataset_list_file.close()


if __name__ == "__main__":
  main()
