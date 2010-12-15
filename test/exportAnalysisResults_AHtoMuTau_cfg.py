import FWCore.ParameterSet.Config as cms

from TauAnalysis.Configuration.recoSampleDefinitionsAHtoMuTau_7TeV_grid_cfi import *

#--------------------------------------------------------------------------------
# Export analysis results into ASCII files,
# in the format used by the CDF collaboration
#--------------------------------------------------------------------------------

process = cms.Process('exportAnalysisResults')

process.DQMStore = cms.Service("DQMStore")

process.maxEvents = cms.untracked.PSet(            
    input = cms.untracked.int32(0)         
)

process.source = cms.Source("EmptySource")

higgsMassPoints = [ '90', '100', '130', '160', '200', '250', '350' ]

process.loadAnalysisResults = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(
        inputFileNames = cms.vstring(
            ##'/data1/veelken/CMSSW_3_8_x/plots/MSSM_Higgs_combined/plotsAHtoMuTau_skimmed.root'
            '/data1/veelken/CMSSW_3_8_x/plots/MSSM_Higgs_combined/plotsAHtoMuTau_TaNCloose_skimmed.root'
        ),
        dqmDirectory_store = cms.string('/export')
    )
)

process.dumpDQMStore = cms.EDAnalyzer("DQMStoreDump")

AHtoMuTau_systematics = [
##    'muonPtUp',
##    'muonPtDown',
##    'tauJetEnUp',
##    'tauJetEnDown',
##    'jetEnUp',
##    'jetEnDown'
]

dqmDirectoryTemplateAHtoMuTau_woBtag = \
  '/export/harvested/%s/ahMuTauAnalyzerOS_woBtag/afterEvtSelDiTauCandidateForAHtoMuTauZeroCharge/'
dqmDirectoryFilterStatAHtoMuTau_woBtag = \
  '/export/harvested/%s/ahMuTauAnalyzerOS_woBtag/FilterStatistics/'
dqmDirectoryTemplateAHtoMuTau_wBtag = \
  '/export/harvested/%s/ahMuTauAnalyzerOS_wBtag/afterEvtSelDiTauCandidateForAHtoMuTauZeroCharge/'
dqmDirectoryFilterStatAHtoMuTau_wBtag = \
  '/export/harvested/%s/ahMuTauAnalyzerOS_wBtag/FilterStatistics/'
meNameTemplateAHtoMuTau = 'DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass'
meNameNumEventsProcessedAHtoMuTau = 'genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1'
meNameNumEventsPassedAHtoMuTau_woBtag = 'evtSelDiTauCandidateForAHtoMuTauZeroCharge/passed_cumulative_numWeighted#a1#s1'
meNameNumEventsPassedAHtoMuTau_wBtag = 'evtSelDiTauCandidateForAHtoMuTauZeroCharge/passed_cumulative_numWeighted#a1#s1'

process.exportAnalysisResults = cms.EDAnalyzer("DQMExportAnalysisResults",

    channels = cms.VPSet(
        cms.PSet(
            name = cms.string("AHtoMuTau_woBtag"),
            shortName = cms.string("ma"),
            binning = cms.string(
                dqmDirectoryTemplateAHtoMuTau_woBtag % 'ZtautauSum'
               + meNameTemplateAHtoMuTau
            ),
            dataIntLumi = cms.double(ZtoMuTau.TARGET_LUMI)
        ),
        cms.PSet(
            name = cms.string("AHtoMuTau_wBtag"),
            shortName = cms.string("mab"),
            binning = cms.string(
                dqmDirectoryTemplateAHtoMuTau_wBtag % 'ZtautauSum'
               + meNameTemplateAHtoMuTau
            ),
            dataIntLumi = cms.double(ZtoMuTau.TARGET_LUMI)
        )
    ),

    outputFilePath = cms.string("/data1/veelken/CMSSW_3_8_x/plots/export_AHtoMuTau_TaNCloose"),                 

    processes = cms.PSet(
        Ztautau = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau_woBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_woBtag % 'ZtautauSum'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'ZtautauSum'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'ZtautauSum'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_woBtag
                         )
                    )
                ),
                AHtoMuTau_wBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_wBtag % 'ZtautauSum'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'ZtautauSum'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'ZtautauSum'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_wBtag
                         )
                    )
                )
            ),
            xSection = cms.double(ZtoMuTau.RECO_SAMPLES['ZtautauPU156bx']['x_sec']
                                 + ZtoMuTau.RECO_SAMPLES['qqZtautauPU156bx']['x_sec']),
            outputFilePath = cms.string(""),
	    outputFileName = cms.string("ztt_#CHANNEL_OUTPUTFILENAME#.hst")
        ),
        Zmumu = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau_woBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_woBtag % 'Zmumu'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'Zmumu'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'Zmumu'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_woBtag
                         )
                    )
                ),
                AHtoMuTau_wBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_wBtag % 'Zmumu'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'Zmumu'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'Zmumu'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_wBtag
                         )
                    )
                )
            ),
            xSection = cms.double(ZtoMuTau.RECO_SAMPLES['ZtautauPU156bx']['x_sec']),
            outputFilePath = cms.string(""),                                   
	    outputFileName = cms.string("zmm_#CHANNEL_OUTPUTFILENAME#.hst")
        ),                                      
        QCD = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau_woBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_woBtag % 'qcdSum'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'qcdSum'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'qcdSum'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_woBtag
                         )
                    )
                ),
                AHtoMuTau_wBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_wBtag % 'qcdSum'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'qcdSum'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'qcdSum'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_wBtag
                         )
                    )
                )
            ),
            outputFilePath = cms.string(""),
	    outputFileName = cms.string("qcd_#CHANNEL_OUTPUTFILENAME#.hst")
        ),    
        WplusJets = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau_woBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_woBtag % 'WplusJets'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'WplusJets'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'WplusJets'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_woBtag
                         )
                    )
                ),
	        AHtoMuTau_wBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_wBtag % 'WplusJets'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'WplusJets'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'WplusJets'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_wBtag
                         )
                    )
                )
            ),
            xSection = cms.double(ZtoMuTau.RECO_SAMPLES['Wenu']['x_sec']
                                 + ZtoMuTau.RECO_SAMPLES['Wmunu']['x_sec']
                                 + ZtoMuTau.RECO_SAMPLES['Wtaunu']['x_sec']),
            outputFilePath = cms.string(""),
	    outputFileName = cms.string("wjets_#CHANNEL_OUTPUTFILENAME#.hst")
        ),
        TTplusJets = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau_woBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_woBtag % 'TTplusJets'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'TTplusJets'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'TTplusJets'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_woBtag
                         )
                    )
                ),
	        AHtoMuTau_wBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_wBtag % 'TTplusJets'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'TTplusJets'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'TTplusJets'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_wBtag
                         )
                    )
                )
            ),
            xSection = cms.double(ZtoMuTau.RECO_SAMPLES['TTplusJets']['x_sec']),
            outputFilePath = cms.string(""),
	    outputFileName = cms.string("ttbar_#CHANNEL_OUTPUTFILENAME#.hst")
        ),
        data = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau_woBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_woBtag % 'data'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'data'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_woBtag % 'data'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_woBtag
                         )
                    )
                ),
		AHtoMuTau_wBtag = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau_wBtag % 'data'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'data'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau_wBtag % 'data'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_wBtag
                         )
                    )
                )
            ),
            outputFilePath = cms.string(""),
	    outputFileName = cms.string("data_#CHANNEL_OUTPUTFILENAME#.hst")
        )
    )

    ##systematics = cms.PSet(
    ##	  muonPtUp = cms.PSet(
    ##	      dqmDirectory = cms.string('sysUncertaintyHistManagerResults/sysMuonPtUp'),
    ## 	      outputFilePath = cms.string("mu_pt+1")
    ##    ),
    ##	  muonPtDown = cms.PSet(
    ##	      dqmDirectory = cms.string('sysUncertaintyHistManagerResults/sysMuonPtDown'),
    ##	      outputFilePath = cms.string("mu_pt-1")
    ##    ),                                                  
    ##	  tauJetEnUp = cms.PSet(
    ##	      dqmDirectory = cms.string('sysUncertaintyHistManagerResults/sysTauJetEnUp'),
    ## 	      outputFilePath = cms.string("tau_es+1")
    ##    ),
    ##	  tauJetEnDown = cms.PSet(
    ##	      dqmDirectory = cms.string('sysUncertaintyHistManagerResults/sysTauJetEnDown'),
    ##	      outputFilePath = cms.string("tau_es-1")
    ##    ),
    ##	  jetEnUp = cms.PSet(
    ##	      dqmDirectory = cms.string('sysUncertaintyHistManagerResults/sysJetEnUp'),
    ## 	      outputFilePath = cms.string("jet_es+1")
    ##    ),
    ##	  jetEnDown = cms.PSet(
    ##	      dqmDirectory = cms.string('sysUncertaintyHistManagerResults/sysJetEnDown'),
    ##	      outputFilePath = cms.string("jet_es-1")
    ##    )                                                
    ##)
)

for higgsMassPoint in higgsMassPoints:

    pset = cms.PSet(
        distributions = cms.PSet(
            AHtoMuTau_woBtag = cms.PSet(
                template = cms.string(
                    dqmDirectoryTemplateAHtoMuTau_woBtag % ('A%sSum' % higgsMassPoint)
                   + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                ),
                systematics = cms.vstring(AHtoMuTau_systematics),
                normalization = cms.PSet(
                    numEventsProcessed = cms.string(
                        dqmDirectoryFilterStatAHtoMuTau_woBtag % ('A%sSum' % higgsMassPoint)
                       + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                    ),
                    numEventsPassed = cms.string(
                        dqmDirectoryFilterStatAHtoMuTau_woBtag % ('A%sSum' % higgsMassPoint)
                       + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_woBtag
                    )
                )
            ),
	    AHtoMuTau_wBtag = cms.PSet(
                template = cms.string(
                    dqmDirectoryTemplateAHtoMuTau_wBtag % ('A%sSum' % higgsMassPoint)
                   + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                ),
                systematics = cms.vstring(AHtoMuTau_systematics),
                normalization = cms.PSet(
                    numEventsProcessed = cms.string(
                        dqmDirectoryFilterStatAHtoMuTau_wBtag % ('A%sSum' % higgsMassPoint)
                       + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                    ),
                    numEventsPassed = cms.string(
                        dqmDirectoryFilterStatAHtoMuTau_wBtag % ('A%sSum' % higgsMassPoint)
                       + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau_wBtag
                    )
                )
            )
        ),
        outputFilePath = cms.string("m%s" % higgsMassPoint),
        outputFileName = cms.string("A_#CHANNEL_OUTPUTFILENAME#.hst")
    )

    xSection_gg = AHtoMuTauSpecific_RECO_SAMPLES['A%s' % higgsMassPoint]['x_sec']
    xSection_bb = AHtoMuTauSpecific_RECO_SAMPLES['bbA%s' % higgsMassPoint]['x_sec']
    setattr(pset, "xSection", cms.double(xSection_gg + xSection_bb))
    
    setattr(process.exportAnalysisResults.processes, "A%s" % higgsMassPoint, pset)

print("computing AHtoMuTau effectice cross-sections...")

processEntryAHtoMuTau = cms.PSet(
    efficiency = cms.PSet(
        numerator = cms.string(meNameNumEventsPassedAHtoMuTau_woBtag),
        denominator = cms.string(meNameNumEventsProcessedAHtoMuTau)
    ),
    numEventsPassed = cms.string(meNameNumEventsPassedAHtoMuTau_woBtag),
    dqmDirectory = cms.string('')
)

process.compDQMEffXsecAHtoMuTau = cms.EDAnalyzer("DQMEffXsecCalculator",
    dataIntLumi = cms.double(ZtoMuTau.TARGET_LUMI),
    channels = cms.PSet(
        QCD = processEntryAHtoMuTau.clone(
            dqmDirectory = cms.string(dqmDirectoryFilterStatAHtoMuTau_woBtag % "qcdSum")
        )
    )
)

process.p = cms.Path(
    process.loadAnalysisResults
   + process.dumpDQMStore 
   + process.exportAnalysisResults
   + process.compDQMEffXsecAHtoMuTau 
)

# print-out all python configuration parameter information
print process.dumpPython()
