import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Phase2C9_cff import Phase2C9
process = cms.Process('HARVESTING',Phase2C9)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D88Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.DQMSaverAtRunEnd_cff')
process.load('Configuration.StandardSequences.Harvesting_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring('file:step3_inDQM.root')
)

# process.options = cms.untracked.PSet(
#     FailPath = cms.untracked.vstring(),
#     IgnoreCompletely = cms.untracked.vstring(),
#     Rethrow = cms.untracked.vstring('ProductNotFound'),
#     SkipEvent = cms.untracked.vstring(),
#     allowUnscheduled = cms.obsolete.untracked.bool,
#     canDeleteEarly = cms.untracked.vstring(),
#     emptyRunLumiMode = cms.obsolete.untracked.string,
#     eventSetup = cms.untracked.PSet(
#         forceNumberOfConcurrentIOVs = cms.untracked.PSet(
#             allowAnyLabel_=cms.required.untracked.uint32
#         ),
#         numberOfConcurrentIOVs = cms.untracked.uint32(1)
#     ),
#     fileMode = cms.untracked.string('FULLMERGE'),
#     forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
#     makeTriggerResults = cms.obsolete.untracked.bool,
#     numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
#     numberOfConcurrentRuns = cms.untracked.uint32(1),
#     numberOfStreams = cms.untracked.uint32(0),
#     numberOfThreads = cms.untracked.uint32(1),
#     printDependencies = cms.untracked.bool(False),
#     sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
#     throwIfIllegalParameter = cms.untracked.bool(True),
#     wantSummary = cms.untracked.bool(False)
# )

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step4 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T15', '')

process.postValidationOuterTracker_step = cms.Path(process.postValidationOuterTracker)
from DQM.SiOuterTracker.OuterTrackerClientConfig_cff import *
process.MyDQMHarvester = cms.Sequence( OuterTrackerClient )
process.DQMHarvestOuterTracker_step = cms.Path(process.MyDQMHarvester)
process.dqmsave_step = cms.Path(process.DQMSaver)

process.schedule = cms.Schedule(process.postValidationOuterTracker_step,process.DQMHarvestOuterTracker_step,process.dqmsave_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
