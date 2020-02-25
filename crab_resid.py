from WMCore.Configuration import Configuration

config = Configuration()
config.section_("General")
config.General.requestName   = 'Layer3_Resid_filter3_JetHT_321833_lowpt_mar20'
config.General.transferLogs = True
config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'python/resid_filter.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.outputFiles = ['PixelTree_1_10.root']
config.section_("Data")
config.Data.inputDataset = '/JetHT/Run2018D-v1/RAW'
config.Data.inputDBS = 'global'
#config.Data.inputDBS = 'phys03'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 1
config.Data.publication = True
config.Data.outLFNDirBase = '/store/user/oamram/PixelTrees'
config.Data.runRange = '321833'
#config.Data.runRange = '324318'
#config.Data.outputDatasetTag  = 'RECO'
#config.Data.outputPrimaryDataset  = 'RelValTTbarLepton_13_RAW'
config.Data.ignoreLocality = True
config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC' 
config.Site.whitelist = ['T1_US_FNAL']
config.Site.ignoreGlobalBlacklist = True
