# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: -s RAW2DIGI,L1Reco,RECO --data --scenario pp --conditions 101X_dataRun2_Express_v7 --era Run2_2018 --eventcontent RECO --datatier RECO --filein /store/express/Run2018A/ExpressPhysics/FEVT/Express-v1/000/315/689/00000/0006526E-0A4F-E811-9FEC-FA163E735BE1.root --python_filename=Data_101X_cfg.py --runUnscheduled -n 10 --no_exec

orbit_begin_array = [0]
orbit_end_array = [1434229500]
index = 0

import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('NTUPLE',eras.Run2_2018)

#--------------------------------------------------------------------------------------------------


### Events to process: 'maxEvents' is already registered by the framework
opt.setDefault('maxEvents', -1)
opt.parseArguments()

# input tracks
# Refitter
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
trackInput = 'generalTracks'
if opt.dataTier == 'RECO' or opt.dataTier == 'FEVT' or opt.dataTier == 'RAW':
    trackInput = 'TrackRefitter'
    process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
    process.TrackRefitter_step = cms.Path(process.MeasurementTrackerEvent*process.TrackRefitter)

#------------------------------------------
#  Define your Analyzer(s) here
#------------------------------------------

# Efficiency
process.TTRHBuilderAngleAndTemplate.PixelCPE = cms.string('PixelCPEClusterRepair') 

# FPix Resolution
process.FPixResolution_Template = cms.EDAnalyzer('Triplets_FPix',
    triggerSource = cms.InputTag('TriggerResults::HLT'),
    ttrhBuilder = cms.string('WithAngleAndTemplate'),
    doBPix = cms.bool(False),
    doFPix = cms.bool(True),
    orbit_beginning = cms.int32(orbit_begin_array[index]),
    orbit_end = cms.int32(orbit_end_array[index]),
)
process.FPixResolution_Generic = process.FPixResolution_Template.clone(
    ttrhBuilder = cms.string('WithTrackAngle'),
    orbit_beginning = cms.int32(orbit_begin_array[index]),
    orbit_end = cms.int32(orbit_end_array[index]),
)

# TFileService used for both BPix/FPix resolution
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('FPix_resolution.root'),
)


# Paths
process.PhaseIPixelNtuplizer_step = cms.Path(process.PhaseINtuplizerPlugin)
process.Efficiency_step = cms.Path(process.PhaseINtuplizerPlugin)

#process.Efficiency_step     = cms.Path(process.TimingStudy)
process.FPixResolution_step = cms.Path(process.FPixResolution_Template)

#------------------------------------------
#  Configurations from cmsDriver.py
#------------------------------------------

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

# number of events to run on
process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(opt.maxEvents)
    input = cms.untracked.int32(20)    
)

# MessageLogger
process.MessageLogger.cerr.FwkReport.reportEvery = 10 if opt.dataTier == 'RAW' else 100

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'/store/express/Run2018A/ExpressPhysics/FEVT/Express-v1/000/315/689/00000/0006526E-0A4F-E811-9FEC-FA163E735BE1.root'
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/321/149/00000/A44B33BD-DD9D-E811-88BC-FA163EF9728B.root'                                  
#'/store/data/Run2018B/SingleMuon/RAW/v1/000/319/077/00000/D69F27EA-887D-E811-BB02-FA163E1FD66F.root'
#'/store/data/Run2017G/ZeroBias/RAW/v1/000/306/521/00000/8434DA8B-ADC7-E711-B9D0-02163E01A20A.root'
#'/store/data/Commissioning2018/ZeroBias/RAW/v1/000/313/142/00000/3A62F915-3B34-E811-9BC9-FA163E6F382C.root'
#'/store/data/Run2017H/SingleMuon/RAW/v1/000/307/082/00000/FCE9D159-6ED2-E711-AC79-02163E01390F.root'
#'/store/data/Run2017H/SingleMuon/RAW/v1/000/307/082/00000/FCFBBF2D-59D2-E711-808E-02163E0145DB.root'
#'/store/data/Commissioning2018/ZeroBias/RAW/v1/000/313/025/00000/64AD95F5-8732-E811-9F3D-FA163E42C0B5.root'
#'/store/data/Run2018D/ZeroBias/RAW/v1/000/324/410/00000/00739499-AFC9-3147-979D-413F734229EF.root'
#'/store/express/Run2018D/ExpressPhysics/FEVT/Express-v1/000/324/318/00000/16238B89-A4D3-3542-A1BA-AC51A3F00DBD.root'
#'/store/data/Run2018D/SingleMuon/RAW/v1/000/324/318/00000/35BD0831-687D-194F-A07C-3E85CF60D2B1.root'
'/store/data/Run2018D/SingleMuon/RAW/v1/000/321/833/00000/9CD43785-D4A9-E811-893C-FA163EF8F660.root'
                                      ),
                            )
#opt.dataTier = 'FEVT'

if opt.secondaryInputFile != '':
    process.source.secondaryFileNames = cms.untracked.vstring(opt.secondaryInputFile)

process.options = cms.untracked.PSet(

)

# Other statements
# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Prompt_v7', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_PromptLike_v7', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_Prompt_Candidate_2018_10_26_20_13_12', '')
#------------------ LOCAL CONDITIONS

#process.GlobalTag.toGet = cms.VPSet(
#    cms.PSet(record = cms.string('SiPixelGainCalibrationOfflineRcd'),
#             tag = cms.string('SiPixelGainCalibration_2018_v3'),
#             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
#             ),
#    cms.PSet(record = cms.string('SiPixelTemplateDBObjectRcd'),
#             tag = cms.string('SiPixelTemplateDBObject_phase1_38T_2018_v6'),
#             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
#             ),
#    cms.PSet(record = cms.string('SiPixelGenErrorDBObjectRcd'),
#             tag = cms.string('SiPixelGenErrorDBObject_phase1_38T_2018_v6'),
#             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
#             ),
#    cms.PSet(record = cms.string('SiPixelLorentzAngleRcd'),
#             tag = cms.string('SiPixelLorentzAngle_phase1_2018_v6'),
#             connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
#             ),
#    cms.PSet(record = cms.string('SiPixel2DTemplateDBObjectRcd'),
#            tag = cms.string('SiPixel2DTemplateDBObject_phase1_38T_2018_v11_num'),
#            connect = cms.string("frontier://FrontierProd/CMS_CONDITIONS"),
#             ),
#    )

#------------------------------------------------------------------------------
# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)

# Schedule definition
if opt.dataTier == 'RECO' or opt.dataTier == 'FEVT':
    process.schedule = cms.Schedule()
else:
    process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step)

# Add TrackRef tter if needed
trackInput == 'TrackRefitter'
if trackInput == 'TrackRefitter': process.schedule.append(process.TrackRefitter_step)
# Add ntuplizers
process.schedule.append(process.FPixResolution_step)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('-s nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition


# Additional output definition
# Schedule definition
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
