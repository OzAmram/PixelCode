# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RECO -s RAW2DIGI,L1Reco,RECO --data --scenario pp --conditions 92X_dataRun2_Prompt_forTier0Replay_PixelLocalReco_TkAl_newVCal_v2 --era Run2_2017 --process NTUPLE --eventcontent RECO --datatier RECO --filein /store/data/Tier0_REPLAY_vocms015/ZeroBias/RECO/PromptReco-v141/00000/6658E50E-9A65-E711-A4F9-02163E011C17.root --python_filename=run_Resolution_ReReco_Data_92X_cfg.py --runUnscheduled -n 10
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

index = 0

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
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Prompt_Candidate_2018_10_26_20_13_12', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.endjob_step,process.RECOoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)









# Insert resolution stuff here

# Refitter
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.TrackRefitter_step = cms.Path(process.MeasurementTrackerEvent*process.TrackRefitter)

#------------------------------------------
#  Define your Analyzer(s) here
#------------------------------------------

# BPix Resolution
process.BPixResolution_Template = cms.EDAnalyzer('Pixel',
    triggerSource = cms.InputTag('TriggerResults::HLT'),
    ttrhBuilder = cms.string('WithAngleAndTemplate'),
    orbit_beginning = cms.int32(orbit_begin_array[index]),
    orbit_end = cms.int32(orbit_end_array[index]),
)
process.BPixResolution_Generic = process.BPixResolution_Template.clone(
    ttrhBuilder = cms.string('WithTrackAngle'),
)

# FPix Resolution
process.FPixResolution_Template = cms.EDAnalyzer('Pixel_phase1',
    triggerSource = cms.InputTag('TriggerResults::HLT'),
    ttrhBuilder = cms.string('WithAngleAndTemplate'),
    doBPix = cms.bool(False),
    doFPix = cms.bool(True)
)
process.FPixResolution_Generic = process.FPixResolution_Template.clone(
    ttrhBuilder = cms.string('WithTrackAngle')
)

# TFileService used for both BPix/FPix resolution
process.TFileService = cms.Service('TFileService',
    fileName = cms.string("BPix_resid_2D_var.root"),
)

# Paths
#process.TTRHBuilderAngleAndTemplate.PixelCPE = cms.string('PixelCPEClusterRepair') 
process.BPixResolution_step = cms.Path(process.BPixResolution_Template)
process.FPixResolution_step = cms.Path(process.FPixResolution_Template*process.FPixResolution_Generic)

# Schedule definition
process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.L1Reco_step,
    process.reconstruction_step,
    process.TrackRefitter_step,
#    process.FPixResolution_step,
    process.BPixResolution_step
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
