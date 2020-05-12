#config file for resid tool. Calculates residuals on layer1 with refitting tracks without layer1 hits
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('NTUPLE',eras.Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

orbit_begin_array = [0]
orbit_end_array = [1434229500]


layer = 1

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
 
#"/store/data/Run2017C/SingleMuon/RECO/PromptReco-v3/000/300/780/00000/044C579F-7E7E-E711-B516-02163E0126FD.root",
#'/store/data/Run2018D/L1Accept/RAW/v1/000/324/420/00000/276128BC-58C0-D14D-9AC4-F2A979FF49E3.root'
#'/store/data/Run2018D/SingleMuon/RAW/v1/000/324/318/00000/35BD0831-687D-194F-A07C-3E85CF60D2B1.root'
#'/store/express/Run2018D/ExpressPhysics/FEVT/Express-v1/000/324/318/00000/16238B89-A4D3-3542-A1BA-AC51A3F00DBD.root'
#'/store/data/Run2018D/SingleMuon/RAW/v1/000/324/318/00000/35BD0831-687D-194F-A07C-3E85CF60D2B1.root'
#'/store/data/Run2018D/JetHT/RAW/v1/000/324/318/00000/53544E6D-FD2B-F843-B87C-43998833DC43.root'
#'/store/data/Run2018D/DoubleMuon/RAW/v1/000/321/833/00000/8644BCBB-D8A9-E811-B640-FA163EDE417A.root'
#'/store/data/Run2018D/DoubleMuon/RAW/v1/000/321/833/00000/B62DB1BB-D9A9-E811-BC1F-FA163E7C3F50.root'
'/store/data/Run2018D/SingleMuon/RAW/v1/000/321/833/00000/9CD43785-D4A9-E811-893C-FA163EF8F660.root'
#'/store/data/Commissioning2018/MinimumBias/RAW/v1/000/312/309/00000/9032B062-F52C-E811-A98B-FA163EBF57CC.root'

),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('RECO nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('RECO_RAW2DIGI_L1Reco_RECO.root'),
    outputCommands = process.RECOEventContent.outputCommands,
    splotLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Prompt_v7', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_PromptLike_v7', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Prompt_Candidate_2018_10_23_14_42_31', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Prompt_Candidate_2018_10_26_20_13_12', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '105X_dataRun2_v6', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)







# Insert resolution stuff here

# Refitter
process.load("RecoTracker.FinalTrackSelectors.TrackerTrackHitFilter_cff")
process.TrackerTrackHitFilter.src = 'generalTracks'
if(layer == 1): process.TrackerTrackHitFilter.commands = cms.vstring("drop PXB","keep PXB 2","keep PXB 3","keep PXB 4","keep PXE","keep TIB","keep TID","keep TOB","keep TEC")
if(layer == 2): process.TrackerTrackHitFilter.commands = cms.vstring("keep PXB","drop PXB 2","keep PXB 3","keep PXB 4","keep PXE","keep TIB","keep TID","keep TOB","keep TEC")
if(layer == 3): process.TrackerTrackHitFilter.commands = cms.vstring("keep PXB","keep PXB 2","drop PXB 3","keep PXB 4","keep PXE","keep TIB","keep TID","keep TOB","keep TEC")

#process.TrackerTrackHitFilter.commands = cms.vstring("keep PXB","keep PXB 2","keep PXB 3","drop PXB 4","keep PXE","keep TIB","keep TID","keep TOB","keep TEC")
#Refit tracks after hit filter
#You might want to customize the options, or use a different refitter
import RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cff
process.trackFitter = RecoTracker.TrackProducer.CTFFinalFitWithMaterial_cff.ctfWithMaterialTracks.clone()
process.trackFitter.src = 'TrackerTrackHitFilter'
#process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.TrackRefitter_step = cms.Path( process.MeasurementTrackerEvent * process.TrackerTrackHitFilter * process.trackFitter * process.reconstruction) 

#TURN ON CR FOR RECO

#process.TTRHBuilderAngleAndTemplate.PixelCPE = cms.string('PixelCPEClusterRepair') 


#------------------------------------------
#  Define your Analyzer(s) here
#------------------------------------------
#process.Layer1_Residuals = cms.EDAnalyzer('Phase1PixelNtuplizer2',
#    trajectoryInput = cms.InputTag('TrackRefitter'))

# BPix Resolution
process.Refit_Residuals = cms.EDAnalyzer('Resid_filter',
        triggerSource = cms.InputTag('TriggerResults::HLT'),
        ttrhBuilder = cms.string('WithAngleAndTemplate'),
        trackInput = cms.InputTag('trackFitter'),
        trackInputGeneral = cms.InputTag('generalTracks'),
        dropLayer = cms.int32(layer)
)

# TFileService used for both BPix/FPix resolution
process.TFileService = cms.Service('TFileService',
    fileName = cms.string("Refit_residuals_test.root"),
)

# Paths
process.Refit_Residuals_step = cms.Path(process.Refit_Residuals)

process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.L1Reco_step,
    process.reconstruction_step,
    process.TrackRefitter_step,
    process.Refit_Residuals_step
    )

# end of insert to cmsDriver script


#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
