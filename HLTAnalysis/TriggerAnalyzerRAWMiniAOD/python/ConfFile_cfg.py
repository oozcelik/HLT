import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
   'file:TEST_HLT.root')
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.demo = cms.EDAnalyzer('TriggerAnalyzerRAWMiniAOD'
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('test.root')
)

process.p = cms.Path(process.demo)
