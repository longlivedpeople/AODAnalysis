import FWCore.ParameterSet.Config as cms


AODanalyzer = cms.EDAnalyzer('AODAnalysis',
    nameOfOutput = cms.string('outputAOD.root'),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    GeneralTrackCollection = cms.InputTag("generalTracks"),
    GenParticleCollection = cms.InputTag("genParticles"),
    PhotonCollection = cms.InputTag("photons"),
    PrimaryVertexCollection = cms.InputTag("offlinePrimaryVertices"),
    bits = cms.InputTag("TriggerResults","","HLT"),
    triggerSummary = cms.InputTag("hltTriggerSummaryAOD")
)


