#!/usr/bin/env python

import sys, os, string, re, copy
from optparse import OptionParser
from ROOT import TFile


# usage description
usage = """Usage: ./mergeDatasets.py [options]\n
Example: ./mergeDatasets.py -w LXBatch_Jobs -d datasetListForMerging.txt\n
For more help: ./mergeDatasets.py  --help
"""

def main():

  # input parameters
  parser = OptionParser(usage=usage)

  parser.add_option("-w", "--main_workdir", dest="main_workdir",
                    help="Main working directory",
                    metavar="MAIN_WORKDIR")

  parser.add_option("-d", "--dataset_list_for_merging", dest="dataset_list_for_merging",
                    help="List of datasets to be merged",
                    metavar="DATASET_LIST_FOR_MERGING")

  parser.add_option("-o", "--output_dir", dest="output_dir",
                    help="Output directory (This parameter is optional)",
                    metavar="OUTPUT_DIR")

  parser.add_option("-a", "--analyzer_module", dest="analyzer_module",
                    help="Name of the analyzer module: SemiLeptanic, Hadronic",
                    default='SemiLeptanic',
                    metavar="ANALYZER_MODULE")

  (options, args) = parser.parse_args()

  # make sure all necessary input parameters are provided
  if not (options.main_workdir and options.dataset_list_for_merging):
    print usage
    sys.exit()

  main_workdir = options.main_workdir
  # redefine main_workdir as an absolute path (if not defined in such form already)
  if not re.search('^/', main_workdir):
    main_workdir = os.path.join(os.getcwd(),main_workdir)

  output_dir = main_workdir
  if options.output_dir:
    output_dir = options.output_dir
    # redefine output_dir as an absolute path (if not defined in such form already)
    if not re.search('^/', output_dir):
      output_dir = os.path.join(os.getcwd(),output_dir)

  group = ''
  group_datasets = {}
  group_xs = {}
  group_L = {}
  dataset_xs = {}

  # open and read the dataset_list_for_merging
  dataset_list_for_merging = open(options.dataset_list_for_merging,'r')
  dataset_list_for_merging_lines = dataset_list_for_merging.readlines()

  for line in dataset_list_for_merging_lines:

    line_elements = line.split()
    if (len(line_elements)==0 or line_elements[0][0]=='#'): continue

    if re.search(':', line):
      group = line_elements[0].rstrip(':')
      group_datasets[group] = []
      group_xs[group] = 0.
      group_L[group] = -1.
      if re.search('L=', line):
        group_L[group] = float(line.split('L=')[-1].strip('\n'))
    else:
      dataset = line_elements[0]
      xs = float(line_elements[1])
      group_datasets[group].append(dataset)
      if xs > 0.:
        group_xs[group] = group_xs[group] + xs
      else:
        group_xs[group] = -1.
      dataset_xs[dataset] = xs

  # final output file
  filename='Final_histograms'
  if (len(options.analyzer_module)>0):
    filename+=str('_'+options.analyzer_module)
  filename+=str('.root')
  output_root_file = TFile( os.path.join(output_dir,filename), 'RECREATE' )

  # write histograms
  groups = group_datasets.keys()
  groups.sort()

  for group in groups:
    print group
    final_histos = {}

    for dataset in group_datasets[group]:
      input_root_file  = os.path.join(main_workdir,dataset.lstrip('/').replace('/','__') + '.root')
      if not os.path.isfile(input_root_file):
        print 'ERROR: File ' + input_root_file + ' not found'
        print 'Aborting'
        sys.exit(1)

      # open input ROOT file
      root_file = TFile(input_root_file)
      htemp = root_file.Get(os.path.join(options.analyzer_module,'h_cutflow'))
      nEventsAll = htemp.GetBinContent(1)
      nEventsStored = htemp.GetBinContent(1)
      scale = 1.
      if group_xs[group] > 0.:
        if group_L[group] > 0.:
          scale = (dataset_xs[dataset]*group_L[group])/nEventsAll
        else:
          scale = dataset_xs[dataset]/(group_xs[group]*nEventsAll)
        print dataset + ' -- Events: %.0f (all), %.0f (stored); relative xs: %.8E; scale: %.8E'%(nEventsAll,nEventsStored,(dataset_xs[dataset]/group_xs[group]),scale)
      else:
        print dataset + ' -- Events: %.0f (all), %.0f (stored); scale: %.8E'%(nEventsAll,nEventsStored,scale)

      # get the number of histograms
      nHistos = root_file.Get(options.analyzer_module).GetListOfKeys().GetEntries()

      # loop over histograms in the input ROOT file
      for h in range(0, nHistos):
        histoName = root_file.Get(options.analyzer_module).GetListOfKeys()[h].GetName()
        #print histoName
        htemp = root_file.Get(os.path.join(options.analyzer_module,histoName))
        if htemp.InheritsFrom('TH1'):

          if histoName not in final_histos.keys():
              final_histos[histoName] = copy.deepcopy(htemp)
              final_histos[histoName].SetName(group + '__' + histoName)
              final_histos[histoName].Scale(scale)
          else:
            final_histos[histoName].Add(htemp, scale) 

    output_root_file.cd()
    histos = final_histos.keys()
    histos.sort()
    print "Writing histograms..."
    for histo in histos:
      #print "Writing histogram: " , final_histos[histo].GetName()
      final_histos[histo].Write()
    print 'Done'

  output_root_file.Close()

  print ''
  print 'Final histograms file: ' + os.path.join(output_dir,filename)
  print ''


if __name__ == '__main__':
  main()
