from WMCore.Configuration import Configuration
config = Configuration()
from CRABClient.UserUtilities import getUsernameFromSiteDB
config.section_('General')
config.General.requestName = 'EbyEPbPb2015_HIMinimumBias2_Check'
config.General.transferOutputs = True
config.General.transferLogs = True
config.section_('JobType')
config.JobType.outputFiles = ['EbyETree_minbias.root']
config.JobType.pyCfgParams = ['noprint']
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '/afs/cern.ch/work/j/jcastle/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/EbyEAnalysis/test/ebyeAnalysis_MinBias_cfg.py'
config.JobType.maxJobRuntimeMin = 1315
config.JobType.inputFiles = ['/afs/cern.ch/work/j/jcastle/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/TrackingCode/HIRun2015Ana/rootfile/PbPb_MB_NTT.root']
config.section_('Data')
config.Data.inputDataset = '/HIMinimumBias2/HIRun2015-PromptReco-v1/AOD'
config.Data.runRange = '262548-263757'
config.Data.unitsPerJob = 20
config.Data.publication = False
config.Data.splitting = 'LumiBased'
config.Data.outLFNDirBase = '/store/user/jcastle/EbyEPbPb2015_HIMinimumBias2_Check'
config.Data.lumiMask = 'cert_262548-263757_PromptReco_HICollisions15_JSON.txt'
config.section_('User')
config.section_('Site')
config.Site.whitelist = ['T2_US_Vanderbilt']
config.Site.storageSite = 'T2_US_Vanderbilt'

if __name__ == '__main__':

   from CRABAPI.RawCommand import crabCommand
   from httplib import HTTPException

   # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
   # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
   config.General.workArea = 'crab_projects'

   def submit(config):
       try:
           crabCommand('submit', config = config)
       except HTTPException, hte:
           print hte.headers

   #############################################################################################
   ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
   #############################################################################################
submit(config)

config.General.requestName = 'EbyEPbPb2015_HIMinimumBias3_Check'
config.Data.inputDataset = '/HIMinimumBias3/HIRun2015-PromptReco-v1/AOD'
config.Data.outLFNDirBase = '/store/user/jcastle/EbyEPbPb2015_HIMinimumBias3_Check'
submit(config)

config.General.requestName = 'EbyEPbPb2015_HIMinimumBias4_Check'
config.Data.inputDataset = '/HIMinimumBias4/HIRun2015-PromptReco-v1/AOD'
config.Data.outLFNDirBase = '/store/user/jcastle/EbyEPbPb2015_HIMinimumBias4_Check'
submit(config)
