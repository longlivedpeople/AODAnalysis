import FWCore.ParameterSet.Config as cms


process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("MyAnalysis.AODAnalysis.IFCAAODAnalysis_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '94X_mc2017_realistic_v11'  # or some other global tag depending on your CMSSW release and sample. 
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
#        'file:/eos/cms/store/mc/RunIIFall17DRPremix/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/AODSIM/MUOTrackFix_94X_mc2017_realistic_MuonTrackFix_01_ext2-v1/110000/00C6F634-B333-E911-A907-FA163E4C57A4.root'
#    )
#)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/f/fernance/private/Long_Lived_Analysis/CMSSW_9_4_4/src/Samples/DisplacedSUSY/AOD/EXO-DisplacedSUSY_squarkToQuarkChi_MSquark_1500_MChi_494_ctau_160mm_TuneCP5_14TeV_pythia8_DRPREMIX_N10000.root'
    )
)





process.p = cms.Path(process.AODanalyzer)


