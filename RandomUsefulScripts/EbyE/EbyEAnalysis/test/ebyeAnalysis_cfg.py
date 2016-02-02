import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os

process = cms.Process("EbyEAnalyzer")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration.EventContent.EventContentHeavyIons_cff')
process.load('HeavyIonsAnalysis.Configuration.collisionEventSelection_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load("RecoHI.HiEvtPlaneAlgos.HiEvtPlane_cfi")
#process.load("RecoHI.HiEvtPlaneAlgos.hievtplaneflatproducer_cfi")
process.load("RecoHI.HiCentralityAlgos.CentralityBin_cfi")
process.load("Configuration.StandardSequences.ReconstructionHeavyIons_cff")
process.load("HeavyIonsAnalysis.Configuration.hfCoincFilter_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('HeavyIonsAnalysis.EbyEAnalysis.ebyeana_cfi')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '75X_dataRun2_v13', '')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/j/jcastle/HIFlowCorr-PromptReco-AOD/9EAF711E-CD9A-E511-8C9A-02163E0142D9.root'
  )
)

##
## Uncomment the following and comment the subsequent toGet for local db file
##

#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#process.CondDBCommon.connect = "sqlite_file:flatparms.db"
#process.PoolDBESSource2 = cms.ESSource("PoolDBESSource",
#                                      process.CondDBCommon,
#                                      toGet = cms.VPSet(cms.PSet(record = cms.string('HeavyIonRPRcd'),
#                                                                 tag = cms.string('EPFlattening_HIRun2011_v5320devx02_offline')
#                                                                 )
#                                                        )
#                                      )

#process.GlobalTag.toGet.extend([
#        cms.PSet(record = cms.string("HeavyIonRPRcd"),
#                 tag = cms.string('HeavyIonRPRcd_PbPb2011_5320_v01_offline'),
#                 connect = cms.untracked.string("frontier://FrontierProd/CMS_COND_PAT_000")
#                 )
#        ])

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("EbyETree.root")
)

process.ebyeana.Subevent_Standard = cms.untracked.bool(True)
#process.ebyeana.effTable_ = cms.string("/afs/cern.ch/work/j/jcastle/CMSSW_7_5_8_patch2/src/HeavyIonsAnalysis/TrackingCode/HIRun2015Ana/rootfile/PbPb_MB_NTT.root") ##-- Local Jobs
process.ebyeana.effTable_ = cms.string("PbPb_MB_NTT.root") ##-- CRAB Jobs

process.ebyeana.countHisto = cms.untracked.bool(False)

process.ebyeana.etaMax_ = cms.untracked.double(2.4)
process.ebyeana.minpt_ = cms.untracked.double(1.)
process.ebyeana.maxpt_ = cms.untracked.double(8.0)

# Q2 trigger selection (later runs)
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.Q2HLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.Q2HLT.HLTPaths = [
#    "HLT_HIQ2Bottom005_Centrality1030_v2",
#    "HLT_HIQ2Top005_Centrality1030_v2",
#    "HLT_HIQ2Bottom005_Centrality3050_v2",
#    "HLT_HIQ2Top005_Centrality3050_v2",
#    "HLT_HIQ2Bottom005_Centrality5070_v2",
#    "HLT_HIQ2Top005_Centrality5070_v2"
    "HLT_HIQ2*"
]


process.p = cms.Path(#process.collisionEventSelection*
                     process.Q2HLT*
                     process.hfCoincFilter3*
                     process.primaryVertexFilter*
                     process.centralityBin*
                     process.ebyeana)
