import FWCore.ParameterSet.Config as cms


AODanalyzer = cms.EDAnalyzer('AODAnalysis',
    nameOfOutput = cms.string('outputAOD.root'),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    GeneralTrackCollection = cms.InputTag("generalTracks"),
    PrimaryVertexCollection = cms.InputTag("offlinePrimaryVertices")
)


