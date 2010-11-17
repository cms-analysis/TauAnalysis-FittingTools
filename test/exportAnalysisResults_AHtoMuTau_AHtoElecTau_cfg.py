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

process.loadAnalysisResults = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(
        inputFileNames = cms.vstring(
            '/data1/veelken/CMSSW_3_8_x/plots/MSSM_Higgs_combined/plotsAHtoMuTau_all.root',
            '/data1/veelken/CMSSW_3_8_x/plots/MSSM_Higgs_combined/plotsZtoElecTau_all.root'
        ),
        dqmDirectory_store = cms.string('')
    )
)

process.dumpDQMStore = cms.EDAnalyzer("DQMStoreDump")

process.sumAHtoElecTau = cms.EDAnalyzer("DQMHistAdder",
    qcdSum = cms.PSet(
        dqmDirectories_input = cms.vstring(
            '/summed/harvested/qcdBCtoESum/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto',
            '/summed/harvested/qcdEMenrichedSum/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto'
        ),
        dqmDirectory_output = cms.string('/summed/harvested/qcdSum/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto')
    ),
    WplusJetsSum = cms.PSet(
        dqmDirectories_input = cms.vstring(
            '/summed/harvested/WtoENu/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto',
            '/summed/harvested/WtoTauNu/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto'
        ),
        dqmDirectory_output = cms.string('/summed/harvested/WplusJetsSum/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto')
    )                                   

)

process.compAHtoElecTauPrediction = cms.EDAnalyzer("DQMHistScaler",
    config = cms.VPSet(
        cms.PSet(
            meName_input = cms.string(
                 'harvested/A100Sum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                + 'evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1'
            ),
            meName_output = cms.string(
                 'harvested/A100Sum/ahElecTauAnalyzer_prediction/FilterStatistics/' \
                + 'evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
            ),
            meType = cms.string("hist"),
            scaleFactor = cms.double(0.613)
        ),
        cms.PSet(   
            meName_input = cms.string(
                 'harvested/A130Sum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                + 'evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1'
            ),
            meName_output = cms.string(
                 'harvested/A130Sum/ahElecTauAnalyzer_prediction/FilterStatistics/' \
                + 'evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
            ),
            meType = cms.string("hist"),
            scaleFactor = cms.double(0.613)
        ),
        cms.PSet( 
            meName_input = cms.string(
                 'harvested/A160Sum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                + 'evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1'
            ),
            meName_output = cms.string(
                 'harvested/A160Sum/ahElecTauAnalyzer_prediction/FilterStatistics/' \
                + 'evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
            ),
            meType = cms.string("hist"),
            scaleFactor = cms.double(0.613)
        ),
        cms.PSet(                                         
            meName_input = cms.string(
                 'harvested/A180Sum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                + 'evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1'
            ),
            meName_output = cms.string(
                 'harvested/A180Sum/ahElecTauAnalyzer_prediction/FilterStatistics/' \
                + 'evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
            ),
            meType = cms.string("hist"),
            scaleFactor = cms.double(0.613)
        ),
        cms.PSet(  
            meName_input = cms.string(
                 'harvested/A200Sum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                + 'evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1'
            ),
            meName_output = cms.string(
                 'harvested/A200Sum/ahElecTauAnalyzer_prediction/FilterStatistics/' \
                + 'evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
            ),
            meType = cms.string("hist"),
            scaleFactor = cms.double(0.613)
        ),
        cms.PSet(  
            meName_input = cms.string(
                 'harvested/A300Sum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                + 'evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1'
            ),
            meName_output = cms.string(
                 'harvested/A300Sum/ahElecTauAnalyzer_prediction/FilterStatistics/' \
                + 'evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
            ),
            meType = cms.string("hist"),
            scaleFactor = cms.double(0.613)
        )
    )
)                                                   

AHtoMuTau_systematics = [
##    'muonPtUp',
##    'muonPtDown',
##    'tauJetEnUp',
##    'tauJetEnDown',
##    'jetEnUp',
##    'jetEnDown'
]

AHtoElecTau_systematics = [
##    'elecEnUp',
##    'elecEnDown',
##    'tauJetEnUp',
##    'tauJetEnDown',
##    'jetEnUp',
##    'jetEnDown'
]

process.exportAnalysisResults = cms.EDAnalyzer("DQMExportAnalysisResults",

    channels = cms.string(
        'AHtoMuTau',
        'AHtoElecTau'
    ),        

    outputFilePath = cms.string("/data1/veelken/CMSSW_3_8_x/plots/MSSM_Higgs_combined/export"),                 

    processes = cms.PSet(
        Ztautau = cms.PSet
            distributions = cms.PSet(
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        '/harvested/A%sSum/ahMuTauAnalyzer_woBtag/afterEvtSelNonCentralJetEt20bTag/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             '/harvested/A%sSum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         ),
                         numEventsPassed = cms.string(
                             '/harvested/A%sSum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         )
                    )
                ),        
                AHtoElecTau = cms.PSet(
                    template = cms.string(
                        '/summed/harvested/Ztautau/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoElecTau_systematics),
                    normalization = cms.PSet(
                        numEventsProcessed = cms.string(
                            '/summed/harvested/Ztautau/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        ),
                        numEventsPassed = cms.string(
                            '/summed/harvested/Ztautau/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        )
                    )
                )
            ),
	    outputFileName = cms.string("ztt_#CHANNEL_OUTPUTFILENAME#.hst")
        ),
        Zee = cms.PSet(
            distributions = cms.PSet(
                AHtoElecTau = cms.PSet(
                    template = cms.string(
                        '/summed/harvested/Zee/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoElecTau_systematics),
                    normalization = cms.PSet(
                        numEventsProcessed = cms.string(
                            '/summed/harvested/Zee/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        ),
                        numEventsPassed = cms.string(
                            '/summed/harvested/Zee/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        )
                    )
                )
            ),
	    outputFileName = cms.string("zee_#CHANNEL_OUTPUTFILENAME#.hst")
        ),
        Zmumu = cms.PSet
            distributions = cms.PSet(
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        '/harvested/Zmumu/ahMuTauAnalyzer_woBtag/afterEvtSelNonCentralJetEt20bTag/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             '/harvested/Zmumu/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         ),
                         numEventsPassed = cms.string(
                             '/harvested/Zmumu/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         )
                    )
                )
            ),
	    outputFileName = cms.string("zmm_#CHANNEL_OUTPUTFILENAME#.hst")
        ),                                      
        QCD = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        '/harvested/qcdSum/ahMuTauAnalyzer_woBtag/afterEvtSelNonCentralJetEt20bTag/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             '/harvested/qcdSum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         ),
                         numEventsPassed = cms.string(
                             '/harvested/qcdSum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         )
                    )
                ),
                AHtoElecTau = cms.PSet(
                    template = cms.string(
                        '/summed/harvested/qcdSum/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoElecTau_systematics),
                    normalization = cms.PSet(
                        numEventsProcessed = cms.string(
                            '/summed/harvested/Ztautau/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        ),
                        numEventsPassed = cms.string(
                            '/summed/harvested/Ztautau/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        )
                    )
                )
            ),
	    outputFileName = cms.string("qcd_#CHANNEL_OUTPUTFILENAME#.hst")
        ),    
        WplusJets = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        '/harvested/WplusJets/ahMuTauAnalyzer_woBtag/afterEvtSelNonCentralJetEt20bTag/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             '/harvested/WplusJets/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         ),
                         numEventsPassed = cms.string(
                             '/harvested/WplusJets/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         )
                    )
                ),
                AHtoElecTau = cms.PSet(
                    template = cms.string(
                        '/summed/harvested/WplusJetsSum/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoElecTau_systematics),
                    normalization = cms.PSet(
                        numEventsProcessed = cms.string(
                            '/summed/harvested/WplusJetsSum/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        ),
                        numEventsPassed = cms.string(
                            '/summed/harvested/WplusJetsSum/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        )
                    )
                )
            ),
	    outputFileName = cms.string("qcd_#CHANNEL_OUTPUTFILENAME#.hst")
        ),
        WplusJets = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        '/harvested/WplusJets/ahMuTauAnalyzer_woBtag/afterEvtSelNonCentralJetEt20bTag/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             '/harvested/WplusJets/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         ),
                         numEventsPassed = cms.string(
                             '/harvested/WplusJets/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         )
                    )
                ),
                AHtoElecTau = cms.PSet(
                    template = cms.string(
                        '/summed/harvested/WplusJetsSum/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoElecTau_systematics),
                    normalization = cms.PSet(
                        numEventsProcessed = cms.string(
                            '/summed/harvested/WplusJetsSum/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        ),
                        numEventsPassed = cms.string(
                            '/summed/harvested/WplusJetsSum/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        )
                    )
                )
            ),
	    outputFileName = cms.string("wjets_#CHANNEL_OUTPUTFILENAME#.hst")
        ),
        TTplusJets = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        '/harvested/TTplusJets/ahMuTauAnalyzer_woBtag/afterEvtSelNonCentralJetEt20bTag/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             '/harvested/TTplusJets/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         ),
                         numEventsPassed = cms.string(
                             '/harvested/TTplusJets/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         )
                    )
                ),
                AHtoElecTau = cms.PSet(
                    template = cms.string(
                        '/summed/harvested/TTplusJetsSum/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoElecTau_systematics),
                    normalization = cms.PSet(
                        numEventsProcessed = cms.string(
                            '/summed/harvested/TTplusJetsSum/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        ),
                        numEventsPassed = cms.string(
                            '/summed/harvested/TTplusJetsSum/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        )
                    )
                )
            ),
	    outputFileName = cms.string("ttbar_#CHANNEL_OUTPUTFILENAME#.hst")
        ),
        gammaPlusJetsSum = cms.PSet(
            distributions = cms.PSet(
                AHtoElecTau = cms.PSet(
                    template = cms.string(
                        '/summed/harvested/TTplusJetsSum/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(AHtoElecTau_systematics),
                    normalization = cms.PSet(
                        numEventsProcessed = cms.string(
                            '/summed/harvested/TTplusJetsSum/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        ),
                        numEventsPassed = cms.string(
                            '/summed/harvested/TTplusJetsSum/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        )
                    )
                )
            ),
	    outputFileName = cms.string("gammajets_#CHANNEL_OUTPUTFILENAME#.hst")
        ),
        data = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        '/harvested/data/ahMuTauAnalyzer_woBtag/afterEvtSelNonCentralJetEt20bTag/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             '/harvested/data/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         ),
                         numEventsPassed = cms.string(
                             '/harvested/data/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                            + '#SYSTEMATICSDIR#/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                         )
                    )
                ),
                AHtoElecTau = cms.PSet(
                    template = cms.string(
                        '/summed/harvested/Data/zElecTauAnalyzer/afterEvtSelElecTauPairZeeHypothesisVeto/' \
                       + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                    ),
                    systematics = cms.vstring(),
                    normalization = cms.PSet(
                        numEventsProcessed = cms.string(
                            '/summed/harvested/Data/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        ),
                        numEventsPassed = cms.string(
                            '/summed/harvested/Data/zElecTauAnalyzer/FilterStatistics/' \
                           + '#SYSTEMATICSDIR#/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                        )
                    )
                )
            ),
	    outputFileName = cms.string("data_#CHANNEL_OUTPUTFILENAME#.hst")
        )
    )

    ##systematics = cms.PSet(
    ##	  elecEnUp = cms.PSet(
    ##	      dqmDirectory = cms.string('sysUncertaintyHistManagerResults/sysElecEnUp'),
    ## 	      outputFilePath = cms.string("mu_pt+1")
    ##    ),
    ##	  elecEnDown = cms.PSet(
    ##	      dqmDirectory = cms.string('sysUncertaintyHistManagerResults/sysElecEnDown'),
    ##	      outputFilePath = cms.string("mu_pt-1")
    ##    ),
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

higgsMassPoints = [ '100', '130', '160', '200', 300' ]
for higgsMassPoint in higgsMassPoints:

    pset = cms.PSet(
        distributions = cms.PSet(
            AHtoMuTau = cms.PSet(
                template = cms.string(
                    'harvested/A%sSum/ahMuTauAnalyzer_woBtag/afterEvtSelNonCentralJetEt20bTag/' \
                   + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                ),
                systematics = cms.vstring(AHtoMuTau_systematics),
                normalization = cms.PSet(
                    numEventsProcessed = cms.string(
                        'harvested/A%sSum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                       + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                    ),
                    numEventsPassed = cms.string(
                        'harvested/A%sSum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                       + '#SYSTEMATICSDIR#/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                    )
                )
            ),
            AHtoElecTau = cms.PSet(
                template = cms.string(
                    'harvested/A%sSum/ahMuTauAnalyzer_woBtag/afterEvtSelNonCentralJetEt20bTag/' \
                   + '#SYSTEMATICSDIR#/DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass' % higgsMassPoint
                ),
                systematics = cms.vstring(AHtoElecTau_systematics),
                normalization = cms.PSet(
                    numEventsProcessed = cms.string(
                        'harvested/A%sSum/ahMuTauAnalyzer_woBtag/FilterStatistics/' \
                       + '#SYSTEMATICSDIR#/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                    ),
                    numEventsPassed = cms.string(
                        'harvested/A%sSum/ahElecTauAnalyzer_prediction/FilterStatistics/' \
                       + '#SYSTEMATICSDIR#/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1' % higgsMassPoint
                    )
                )
            )
        ),
        outputFilePath = cms.string("m%s" % higgsMassPoint),
        outputFileName = cms.string("A_#CHANNEL_OUTPUTFILENAME#.hst")
    )
    
    setattr(process.exportAnalysisResults.processes, "A%s" % higgsMassPoint, pset)

process.compDQMEffXsec = cms.EDAnalyzer("DQMEffXsecCalculator",
    dataIntLumi = cms.double(ZtoMuTau.TARGET_LUMI),
    channels = cms.PSet(
        AHtoElecTau_QCD = cms.PSet(
            efficiency = cms.PSet(
                numerator = cms.string(
                    'FilterStatistics/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
                ),
                denominator = cms.string(
                    'FilterStatistics/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1'
                )
            ),
            numEventsPassed = cms.string(
                'FilterStatistics/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
            ),
            dqmDirectory = cms.string('harvested/qcdSum/zElecTauAnalyzer')
        ),
        AHtoElecTau_WplusJets = cms.PSet(
            efficiency = cms.PSet(
                numerator = cms.string(
                    'FilterStatistics/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
                ),
                denominator = cms.string(
                    'FilterStatistics/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1'
                )
            ),
            numEventsPassed = cms.string(
                'FilterStatistics/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
            ),
            dqmDirectory = cms.string('harvested/WplusJetsSum/zElecTauAnalyzer')
        ),
        AHtoElecTau_gammaPlusJets = cms.PSet(
            efficiency = cms.PSet(
                numerator = cms.string(
                    'FilterStatistics/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
                ),
                denominator = cms.string(
                    'FilterStatistics/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1'
                )
            ),
            numEventsPassed = cms.string(
                'FilterStatistics/evtSelElecTauPairZeeHypothesisVeto/passed_cumulative_numWeighted#a1#s1'
            ),
            dqmDirectory = cms.string('harvested/gammaPlusSum/zElecTauAnalyzer')
        ),
        AHtoMuTau_QCD = cms.PSet(
            efficiency = cms.PSet(
                numerator = cms.string(
                    'FilterStatistics/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1'
                ),
                denominator = cms.string(
                    'FilterStatistics/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1'
                )
            ),
            numEventsPassed = cms.string(
                'FilterStatistics/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1'
            ),
            dqmDirectory = cms.string('harvested/qcdSum/ahMuTauAnalyzer_woBtag')
        )
    )
)

process.p = cms.Path(
    process.loadAnalysisResults
  #+ process.dumpDQMStore 
   + process.exportAnalysisResults_woBtag
   + process.compDQMEffXsec
)
