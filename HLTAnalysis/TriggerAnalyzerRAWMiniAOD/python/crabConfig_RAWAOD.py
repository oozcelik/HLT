#from WMCore.Configuration import Configuration
#config = Configuration()
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()
config.General.transferOutputs = True
config.JobType.psetName = 'TEST_HLT.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['RAWAOD_out.root']
config.Data.inputDataset = '/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/RunIISummer17DRStdmix-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/AODSIM'
config.Data.secondaryInputDataset = '/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/RunIISummer17DRStdmix-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/GEN-SIM-RAW'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.totalUnits = 2
#config.Data.unitsPerJob = 10
config.Data.outLFNDirBase = '/store/user/oozcelik/CRAB_AOUTPUT/HLTEffNtuples/RunIISummer17DRMC/BdToJpsiKstar'
config.section_('User')
config.section_('Site')
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T2_CH_CERN'
