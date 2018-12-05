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
import FWCore.ParameterSet.VarParsing as opts
opt = opts.VarParsing ('analysis')
opt.register('dataTier',           'RAW',
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.string,
	     'Input data tier: RAW, RECO, FEVT or AOD')

opt.register('inputFile',          '',
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.string,
	     'input file name')

opt.register('secondaryInputFile', '',
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.string,
	     'input file name')

opt.register('Efficiency',         False,
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.bool,
	     'Specify if you want to create Efficiency (TimingStudy) ntuples')

opt.register('LATrees',            False,
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.bool,
	     'Specify if you want to create Lorentz Angle (SiPixelLorentzAngle) ntuples')

opt.register('BPixResolution',     False,
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.bool,
	     'Specify if you want to create BPix Resolution (Pxl) ntuples')

opt.register('FPixResolution',     True,
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.bool,
	     'Specify if you want to create FPix Resolution (Pixel) ntuples')


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

process.PhaseINtuplizerPlugin = cms.EDAnalyzer("PhaseIPixelNtuplizer",
     trajectoryInput = cms.InputTag('TrackRefitter'),
      outputFileName = cms.untracked.string('/eos/cms/store/group/dpg_tracker_pixel/comm_pixel/Monitoring/2018/Efficiency/321149_Tv6Gv3/Efficiency_321149_277.root'),
      # Do not save everything and downscale clusters
      clusterSaveDownscaleFactor     = cms.untracked.int32(10000),
      #eventSaveDownscaleFactor       = cms.untracked.int32(opt.prescale),
      eventSaveDownscaleFactor       = cms.untracked.int32(50),
      saveDigiTree                   = cms.untracked.bool(False),
      aveTrackTree                  = cms.untracked.bool(False),
      saveNonPropagatedExtraTrajTree = cms.untracked.bool(False),  
)

# Lorentz Angle
process.SiPixelLorentzAngle = cms.EDAnalyzer("SiPixelLorentzAngle",
        src = cms.InputTag(trackInput),
        fileName = cms.string('/eos/cms/store/group/dpg_tracker_pixel/comm_pixel/Monitoring/2018/LA/321149_Tv6Gv3/LA_321149_277.root'),
        fileNameFit     = cms.string(""),
        binsDepth       = cms.int32(50),
        binsDrift =     cms.int32(200),
        ptMin = cms.double(3),
        ptMinFPix = cms.double(0.1),
        #in case of MC set this to true to save the simhits (does not work currently, Mixing Module needs to be included correctly)
        simData = cms.bool(False),
        normChi2Max = cms.double(2),
        clustSizeYMin = cms.int32(2),
        residualMax = cms.double(0.01),
        clustChargeMax = cms.double(120000)
)

# BPix Resolution
process.BPixResolution_Template = cms.EDAnalyzer('Pxl',
        triggerSource = cms.InputTag('TriggerResults::HLT'),
        ttrhBuilder = cms.string('WithAngleAndTemplate'),
        orbit_beginning = cms.int32(orbit_begin_array[index]),
        orbit_end = cms.int32(orbit_end_array[index]),
        )

process.BPixResolution_Generic = process.BPixResolution_Template.clone(
    ttrhBuilder = cms.string('WithTrackAngle'),
    orbit_beginning = cms.int32(orbit_begin_array[index]),
    orbit_end = cms.int32(orbit_end_array[index]),
)

# FPix Resolution
process.FPixResolution_Template = cms.EDAnalyzer('Pixel_phase1',
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
process.LATrees_step        = cms.Path(process.SiPixelLorentzAngle)
process.BPixResolution_step = cms.Path(process.BPixResolution_Template*process.BPixResolution_Generic)
process.FPixResolution_step = cms.Path(process.FPixResolution_Template*process.FPixResolution_Generic)

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
'/store/data/Run2018D/SingleMuon/RAW/v1/000/324/318/00000/35BD0831-687D-194F-A07C-3E85CF60D2B1.root'
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
process.GlobalTag = GlobalTag(process.GlobalTag, '103X_dataRun2_PromptLike_v7', '')

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
if opt.Efficiency:     process.schedule.append(process.Efficiency_step)
if opt.LATrees:        process.schedule.append(process.LATrees_step)
if opt.BPixResolution: process.schedule.append(process.BPixResolution_step)
if opt.FPixResolution: process.schedule.append(process.FPixResolution_step)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('-s nevts:10'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

#process.RECOoutput = cms.OutputModule("PoolOutputModule",
#    dataset = cms.untracked.PSet(
#        dataTier = cms.untracked.string('RECO'),
#        filterName = cms.untracked.string('')
#    ),
#    fileName = cms.untracked.string('-s_RAW2DIGI_L1Reco_RECO.root'),
#    outputCommands = process.RECOEventContent.outputCommands,
#    splitLevel = cms.untracked.int32(0)
#)

# Additional output definition

# Other statements
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Express_v7', '')

# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
##process.L1Reco_step = cms.Path(process.L1Reco)
#process.reconstruction_step = cms.Path(process.reconstruction)
##process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RECOoutput_step = cms.EndPath(process.RECOoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,p#rocess.endjob_step,process.RECOoutput_step)
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
