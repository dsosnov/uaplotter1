# https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

import os
fullPath=os.getcwd()
cmsswSrcPath = os.environ['CMSSW_BASE']+'/src'

from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.requestName = 'test'
config.General.transferOutputs = True
config.General.transferLogs = True
config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'HI2016_crab.py'
config.JobType.inputFiles = ['crab_script.sh','crab_start.cxx','crab_makeHistoss.cxx']
config.JobType.inputFiles += [cmsswSrcPath+'/UATree/UADataFormat/lib/libUADataFormat.so',cmsswSrcPath+'/UATree/UADataFormat/eventdict_rdict.pcm']
config.JobType.inputFiles += [cmsswSrcPath+'/TOTEMdataFormat/lib/libTOTEMdataFormat.so', cmsswSrcPath+'/TOTEMdataFormat/src/eventdictT_rdict.pcm']
config.JobType.inputFiles += [cmsswSrcPath+'/uaplotter1/lib/libuaplotter.so',            cmsswSrcPath+'/uaplotter1/src/uaplot_dict_rdict.pcm']
config.JobType.scriptExe = 'crab_script.sh'
config.JobType.outputFiles = ['uaplot_output_histos_noise_PARun2016_Pbp.root']
config.section_('Data')
config.Data.inputDataset = '' # dataset
config.Data.unitsPerJob = 50000
config.Data.splitting = 'EventAwareLumiBased'
config.Data.publication = False
config.Data.runRange = 1
config.section_('User')
config.section_('Site')
config.Site.blacklist = ['T2_UA_KIPT'] #['T0_CH_CERN']
#config.Site.whitelist = ['T2_RU_JINR']
config.Site.storageSite = 'T2_RU_JINR'
