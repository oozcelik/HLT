import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    '/store/mc/RunIISummer17DRStdmix/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/GEN-SIM-RAW/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/16BB5AF9-EAAB-E711-B47F-FA163E196D2C.root'
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('test.root')
)

process.demo = cms.EDAnalyzer('TriggerAnalyzerRAWMiniAOD'
)


process.p = cms.Path(process.demo)
process.schedule = cms.Schedule(process.p)
