import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C11M9_cff import Phase2C11M9

process = cms.Process('RERECO', Phase2C11M9)

# Import standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D98Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.Validation_cff')
process.load('DQMServices.Core.DQMStoreNonLegacy_cff')
process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        "/store/relval/CMSSW_14_0_0_pre2/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU_133X_mcRun4_realistic_v1_STD_2026D98_PU200_RV229-v1/2580000/0b2b0b0b-f312-48a8-9d46-ccbadc69bbfd.root",
        #"/store/relval/CMSSW_14_1_0_pre4/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/140X_mcRun4_realistic_v4_STD_2026D110_noPU-v1/2590000/2556e4d4-97da-4708-9cf2-4e2eecbed0ab.root",
        ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(1),  # Single thread for debugging
    numberOfStreams = cms.untracked.uint32(1),  # Single stream for debugging
    wantSummary = cms.untracked.bool(True)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.MessageLogger.cerr.threshold = 'DEBUG'
process.MessageLogger.cerr.default.limit = 10
process.MessageLogger.debugModules = cms.untracked.vstring('*')

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition
process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('step3_pre4_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Other statements
process.mix.playback = True
process.mix.digitizers = cms.PSet()
for a in process.aliases: delattr(process, a)
process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T21', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Phase2 OT rechit step
process.load('RecoLocalTracker.Phase2TrackerRecHits.Phase2TrackerRecHits_cfi')
process.load('RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi')
process.rechits_step = cms.Path(process.siPhase2RecHits * process.siPixelRecHits)

# DQM modules
process.load('DQM.SiTrackerPhase2.Phase2TrackerDQMFirstStep_cff')
process.load('DQM.SiTrackerPhase2.Phase2OTMonitorRecHit_cfi')
process.otdqm_seq = cms.Sequence(process.trackerphase2DQMSource * process.Phase2OTMonitorRecHit)

process.load('Validation.SiTrackerPhase2V.Phase2TrackerValidationFirstStep_cff')
process.load('Validation.SiTrackerPhase2V.Phase2OTValidateRecHit_cfi')
process.otvalid_seq = cms.Sequence(process.trackerphase2ValidationSource * process.Phase2OTValidateRecHit)

process.dqm_step = cms.Path(process.otdqm_seq * process.stubValidOT)
process.validation_step = cms.Path(process.otvalid_seq)

# Schedule definition
process.schedule = cms.Schedule(process.rechits_step,
                                process.dqm_step,
                                process.validation_step,
                                process.DQMoutput_step)

# Customisation of the process
from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn
process = setCrossingFrameOn(process)
