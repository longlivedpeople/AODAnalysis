""" File to especify the samples that are going to be processed 

>> The structure should appear as follows:
We have a variable dict: samples = {}
samples['process'] = {'name' : 'name',
                      'files':['root file 1',
                               'root file 2',
                                ...] }

"""


#samples['weight0'] = {'name' : 'weight0',
#                             'files' : ['/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/AODAnalysis/outputAOD_noWeight.root']}


#samples['weight0p5'] = {'name' : 'weight0p5',
#                             'files' : ['/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/AODAnalysis/outputAOD_w0p5.root']}


#samples['weight0p7'] = {'name' : 'weight0p7',
#                             'files' : ['/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/AODAnalysis/outputAOD_w0p7.root']}

samples['DYJetsToLLM50'] = {'name' : 'DY M-50',
                            'files' : ['/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/AODAnalysis/outputAODDY.root']}




samples['Signal'] = {'name' : 'Signal',
                     'files' : ['/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/MyAnalysis/AODAnalysis/outputAODSignal.root']}


