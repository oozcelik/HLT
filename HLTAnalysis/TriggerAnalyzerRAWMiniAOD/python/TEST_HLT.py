# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('MYHLT',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('HLTrigger.Configuration.HLT_Test_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/mc/RunIISummer17DRStdmix/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/28D812CB-B7AD-E711-915A-008CFAE44F30.root'
'/store/mc/RunIISummer17DRStdmix/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/82747DA9-D1A9-E711-B840-008CFAC93C08.root'
),
    secondaryFileNames = cms.untracked.vstring(
#'/store/mc/RunIISummer17DRStdmix/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/1C31836F-8FAB-E711-9F64-008CFA111174.root',
#'/store/mc/RunIISummer17DRStdmix/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/C83EA7E8-A3A9-E711-AE3C-008CFAE45058.root',
#'/store/mc/RunIISummer17DRStdmix/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/F07521ED-ACA9-E711-A9F4-008CFA197B54.root'
'/store/mc/RunIISummer17DRStdmix/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/54BF46AD-9CA9-E711-A0DC-008CFAE45170.root',
'/store/mc/RunIISummer17DRStdmix/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/7C1FF52E-ADA9-E711-B9A5-008CFAE45188.root',
'/store/mc/RunIISummer17DRStdmix/BdToJpsiKstar_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/00000/CE980620-ADA9-E711-8E62-008CFAE45308.root'
    )
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('TEST nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('TEST_HLT.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_TSG_2017_12_19_13_49_40', '')

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

process.demo = cms.EDAnalyzer('TriggerAnalyzerRAWMiniAOD')
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( "RAWAOD_out.root" )
                                   )
process.demo_step = cms.EndPath(process.demo)


# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.demo_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions

# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
