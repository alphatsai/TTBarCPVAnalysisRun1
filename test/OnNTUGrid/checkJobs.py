#!/usr/bin/env python

import sys, os, shutil, re, subprocess, commands, time
from optparse import OptionParser

# usage description
usage = """Usage: ./createAndSubmitJobs.py [options]\n
Example: ./createAndSubmitJobs.py -w LXBatch_Jobs -d datasetList.txt -c btagvalidation_cfg.py\n
For more help: ./createAndSubmitJobs.py --help
"""

def main():
  # input parameters
  parser = OptionParser(usage=usage)

  parser.add_option("-w", "--main_workdir", dest="main_workdir", action='store', help="Main working directory", metavar="MAIN_WORKDIR")
  parser.add_option("-o", "--output_filename", dest="output_filename", action='store', default='results', help="Output ROOT filename (Default set to results)", metavar="OUTPUT_FILENAME")
  parser.add_option("-S", "--save_path", dest="save_path", action='store', default='', help="Save files path to copy output files to (This parameter is optional)", metavar="SAVE_PATH")
  parser.add_option("-r", "--resubmit", dest="reSubmit", action='store', default=False, help="Resubmit jobs", metavar="RESUBMIT")

  (options, args) = parser.parse_args()

  # make sure all necessary input parameters are provided
  if not (options.main_workdir ):
    print usage
    sys.exit()

  main_workdir = os.path.join(os.getcwd(), options.main_workdir)
  save_path = options.save_path
  output_filename = options.output_filename

  # create the main working directory
  if not os.path.exists(main_workdir) :
	sys.exit() 

  originPath = os.getcwd()
  os.chdir(main_workdir)
  
  if not os.path.exists("datasetList.txt") :
	sys.exit()

  samples = commands.getoutput("cat datasetList.txt | grep -v '#' | awk '{print $1}' | sed 's/^\///g' | sed 's/\//__/g'").split('\n')
 
  for sample in samples:
	print "==============================================================================================="
	print sample
	jobNum = commands.getoutput("ls -l "+sample+"/input | grep job | wc -l ")	
	logNum = commands.getoutput("ls -l "+sample+"/output | grep .log | wc -l ")
	if save_path != '':
		doneJobs = commands.getoutput("ls -l "+save_path+"/"+sample+" | grep root | awk '{print $9}' | sed 's/"+output_filename+"_\(.*\)\.root/\\1/g'").split('\n')
	else:
		doneJobs = commands.getoutput("ls -l "+sample+"/output"+" | grep root | awk '{print $9}' | sed 's/"+output_filename+"_\(.*\)\.root/\\1/g'").split('\n')

	print "Num Log Files "+logNum+"/"+jobNum
	if len(doneJobs) == 1 and doneJobs[0] == '':
		doneNum = 0
	else:	
		doneNum = len(doneJobs) 
	print "Status(root): "+str(doneNum)+"/"+jobNum

	i = 0
	notDone=[]
	while( i<int(jobNum) ):
		if str(i) not in doneJobs:
			notDone.append(i)
		i += 1
	if len(notDone) != 0 :	
		print "Not done: "
		print notDone
		if options.reSubmit:
			for num in notDone:
				commands.getoutput("mv "+sample+"/output/job_"+str(num)+".log "+sample)
				jobName=sample+"."+str(num)
				os.chdir(sample+"/input/job_"+str(num))
         			print 'Submiting '+jobName
 			        subprocess.call("qsub -q cms "+jobName, shell=True)
 			        #print("qsub -q cms "+jobName)
				#print os.getcwd()
				os.chdir(main_workdir)
	else:
		print "Done!"

  os.chdir(originPath)

#  
#   
#  # copy the dataset list file to the main_workdir
#  shutil.copyfile(dataset_list,os.path.join(main_workdir,'datasetList.txt'))
#
#  # copy the CMSSW cfg file to the cfg_files_dir
#  shutil.copyfile(cmssw_cfg,os.path.join(main_workdir,'CMSSW_cfg.py'))
#
#  # look for pileup distribution files and copy them into main_workdir
#  cfg_dirname = os.path.dirname(cmssw_cfg)
#  if cfg_dirname=='':
#    cfg_dirname = os.getcwd()
#  for filename in os.listdir(cfg_dirname):
#    if not os.path.isfile(os.path.join(cfg_dirname,filename)):
#      continue
#
#  # open and read the dataset_list file
#  dataset_list_file = open(dataset_list,"r")
#  dataset_list_lines = dataset_list_file.readlines()
#
#  # loop over datasets
#  for line in dataset_list_lines:
#    line_elements = line.split()
#    if (len(line_elements)==0 or line_elements[0][0]=='#'): continue
#
#    output_filename = options.output_filename
#    cfg_parameters = ''
#    if( len(line_elements)>3 ):
#      cfg_parameters = line_elements[3]
#      if( line_elements[3].split('=')[0]=='outFilename' ): output_filename = line_elements[3].split('=')[1].replace('.root','')
#      for par in range(4,len(line_elements)):
#        cfg_parameters = cfg_parameters + ' ' + line_elements[par]
#        if( line_elements[par].split('=')[0]=='outFilename' ): output_filename = line_elements[par].split('=')[1].replace('.root','')
#
#    dataset = line_elements[0].lstrip('/').replace('/','__')
#    print 'Processing ' + line_elements[0]
#
#    dataset_workdir = os.path.join(main_workdir,dataset)
#
#    # create the dataset working directory
#    os.mkdir(dataset_workdir)
#    os.mkdir(os.path.join(dataset_workdir,'input'))
#    os.mkdir(os.path.join(dataset_workdir,'output'))
#
#    save_path = ' '
#    if ( options.save_path ):
#      save_path = ''.join([options.save_path,"/",dataset])
#      if not os.path.exists(save_path):
#       commands.getoutput('mkdir -p '+save_path)
#      #proc = subprocess.Popen( [ 'mkdir', '-p', save_path ], stdout = subprocess.PIPE, stderr = subprocess.STDOUT )
#
#    filelist = []
#    process_input_dir(line_elements[2], options.match, filelist)
#
#    ##################################################
#    njobs = line_elements[1]
#    numfiles = int(len(filelist)*float(options.fraction))
#    ijobmax=int(njobs)
#    if ijobmax > numfiles:
#       ijobmax=numfiles
#       print '  Number of jobs requested exceeds the total number of input files.\n  The number of jobs set to ' + str(ijobmax) + ' to match the number of input files'
#    filesperjob = int(numfiles/ijobmax)
#    if numfiles%ijobmax!=0:
#       filesperjob = filesperjob + 1
#       ijobmax = int(numfiles/filesperjob)
#       if numfiles%filesperjob!=0:
#           ijobmax = ijobmax + 1
#       if ijobmax != int(njobs):
#           print '  Could not create the exact number of jobs requested.\n  For optimal job splitting, the number of jobs set to ' + str(ijobmax)
#    #################################################
#    originPath = os.getcwd()
#    ifile = 0
#    for ijob in range(ijobmax):
#       # prepare the list of input files
#       input_files = '    \'file:' + filelist[ifile] + '\''
#       for i in range(filesperjob-1):
#           if ifile>(numfiles-2):
#               break
#           ifile = ifile + 1
#           input_files += (',\n    \'file:' + filelist[ifile] + '\'')
#       ifile = ifile + 1
#
#       jobName = dataset_workdir.split('/')[len(dataset_workdir.split('/'))-1]+'.'+str(ijob)
#       ## create input cfi file
#       os.mkdir(os.path.join(dataset_workdir,'input','job_'+str(ijob)))
#       input_files_cfi = open(os.path.join(dataset_workdir,'input','job_'+str(ijob),'inputFiles_cfi.py'),'w')
#       input_files_cfi.write(re.sub('INPUT_FILES',input_files,cfi_template))
#       input_files_cfi.close()
#
#       ## create Bash script
#       #bash_script = open(os.path.join(dataset_workdir,'input','job_'+str(ijob), 'job.sh'),'w')
#       bash_script = open(os.path.join(dataset_workdir,'input','job_'+str(ijob), jobName),'w')
#       bash_script_content = re.sub('MAIN_WORKDIR',main_workdir,bash_template)
#       bash_script_content = re.sub('DATASET_WORKDIR',dataset_workdir,bash_script_content)
#       bash_script_content = re.sub('JOB_NUMBER',str(ijob),bash_script_content)
#       bash_script_content = re.sub('CFG_PARAMETERS',cfg_parameters,bash_script_content)
#       bash_script_content = re.sub('OUTPUT_FILENAME',output_filename,bash_script_content)
#       bash_script_content = re.sub('SAVE_PATH',save_path,bash_script_content)
#       bash_script.write(bash_script_content)
#       bash_script.close()
#
#       if(not options.no_submission):
#         os.chdir(os.path.join(dataset_workdir,'input','job_'+str(ijob)))
#         #jobName = dataset_workdir.split('/')[len(dataset_workdir.split('/'))-1]+'_'+str(ijob)
#         #os.system('chmod +x job.sh')
#         print 'Submiting '+jobName
#         subprocess.call("qsub -q cms "+jobName, shell=True)
#         os.chdir(originPath)
#         time.sleep(1)
#
#  # close all open files
#  dataset_list_file.close()
#
#
if __name__ == "__main__":
  main()
