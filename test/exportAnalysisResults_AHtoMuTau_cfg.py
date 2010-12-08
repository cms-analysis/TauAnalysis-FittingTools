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

higgsMassPoints = [ '100', '130', '160', '200', '300' ]

process.loadAnalysisResults = cms.EDAnalyzer("DQMFileLoader",
    all = cms.PSet(
        inputFileNames = cms.vstring(
            '/data1/veelken/CMSSW_3_8_x/plots/MSSM_Higgs_combined/plotsAHtoMuTau_skimmed.root'
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

dqmDirectoryTemplateAHtoMuTau = '/export/harvested/%s/ahMuTauAnalyzer_woBtag/afterEvtSelNonCentralJetEt20bTag/'
dqmDirectoryFilterStatAHtoMuTau = '/export/harvested/%s/ahMuTauAnalyzer_woBtag/FilterStatistics/'
meNameTemplateAHtoMuTau = 'DiTauCandidateSVfitQuantities/psKine_MEt_ptBalance/Mass'
meNameNumEventsProcessedAHtoMuTau = 'genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1'
meNameNumEventsPassedAHtoMuTau = 'evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1'

process.exportAnalysisResults = cms.EDAnalyzer("DQMExportAnalysisResults",

    channels = cms.VPSet(
        cms.PSet(
            name = cms.string("AHtoMuTau"),
            shortName = cms.string("ma"),
            binning = cms.string(
                dqmDirectoryTemplateAHtoMuTau % 'Ztautau'
               + meNameTemplateAHtoMuTau
            ),
            dataIntLumi = cms.double(35.0)
        )
    ),

    outputFilePath = cms.string("/data1/veelken/CMSSW_3_8_x/plots/export_AHtoMuTau"),                 

    processes = cms.PSet(
        Ztautau = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau % 'Ztautau'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau % 'Ztautau'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau % 'Ztautau'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau
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
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau % 'Zmumu'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau % 'Zmumu'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau % 'Zmumu'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau
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
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau % 'qcdSum'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau % 'qcdSum'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau % 'qcdSum'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau
                         )
                    )
                )
            ),
            outputFilePath = cms.string(""),
	    outputFileName = cms.string("qcd_#CHANNEL_OUTPUTFILENAME#.hst")
        ),    
        WplusJets = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau % 'WplusJets'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau % 'WplusJets'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau % 'WplusJets'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau
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
        ##TTplusJets = cms.PSet(
        ##    distributions = cms.PSet(
        ##        AHtoMuTau = cms.PSet(
        ##            template = cms.string(
        ##                dqmDirectoryTemplateAHtoMuTau % 'TTplusJets'
        ##               + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
        ##            ),
        ##            systematics = cms.vstring(AHtoMuTau_systematics),
        ##            normalization = cms.PSet(
        ##                 numEventsProcessed = cms.string(
        ##                     dqmDirectoryFilterStatAHtoMuTau % 'TTplusJets'
        ##                    + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
        ##                 ),
        ##                 numEventsPassed = cms.string(
        ##                     dqmDirectoryFilterStatAHtoMuTau % 'TTplusJets'
        ##                    + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau
        ##                 )
        ##            )
        ##        )
        ##    ),
        ##    xSection = cms.double(ZtoMuTau.RECO_SAMPLES['TTplusJets']['x_sec']),
        ##    outputFilePath = cms.string(""),
	##    outputFileName = cms.string("ttbar_#CHANNEL_OUTPUTFILENAME#.hst")
        ##),
        data = cms.PSet(
            distributions = cms.PSet(
                AHtoMuTau = cms.PSet(
                    template = cms.string(
                        dqmDirectoryTemplateAHtoMuTau % 'data'
                       + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                    ),
                    systematics = cms.vstring(AHtoMuTau_systematics),
                    normalization = cms.PSet(
                         numEventsProcessed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau % 'data'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                         ),
                         numEventsPassed = cms.string(
                             dqmDirectoryFilterStatAHtoMuTau % 'data'
                            + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau
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
            AHtoMuTau = cms.PSet(
                template = cms.string(
                    dqmDirectoryTemplateAHtoMuTau % ('A%sSum' % higgsMassPoint)
                   + '#SYSTEMATICSDIR#/' + meNameTemplateAHtoMuTau
                ),
                systematics = cms.vstring(AHtoMuTau_systematics),
                normalization = cms.PSet(
                    numEventsProcessed = cms.string(
                        dqmDirectoryFilterStatAHtoMuTau % ('A%sSum' % higgsMassPoint)
                       + '#SYSTEMATICSDIR#/' + meNameNumEventsProcessedAHtoMuTau
                    ),
                    numEventsPassed = cms.string(
                        dqmDirectoryFilterStatAHtoMuTau % ('A%sSum' % higgsMassPoint)
                       + '#SYSTEMATICSDIR#/' + meNameNumEventsPassedAHtoMuTau
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

numeratorAHtoMuTau = 'FilterStatistics/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1'
denominatorAHtoMuTau = 'FilterStatistics/genPhaseSpaceCut/passed_cumulative_numWeighted#a1#s1'
numEventsPassedAHtoMuTau = 'FilterStatistics/evtSelNonCentralJetEt20bTag/passed_cumulative_numWeighted#a1#s1'
processEntryAHtoMuTau = cms.PSet(
    efficiency = cms.PSet(
        numerator = cms.string(numeratorAHtoMuTau),
        denominator = cms.string(denominatorAHtoMuTau)
    ),
    numEventsPassed = cms.string(numEventsPassedAHtoMuTau),
    dqmDirectory = cms.string('')
)

process.compDQMEffXsecAHtoMuTau = cms.EDAnalyzer("DQMEffXsecCalculator",
    dataIntLumi = cms.double(ZtoMuTau.TARGET_LUMI),
    channels = cms.PSet(
        QCD = processEntryAHtoMuTau.clone(
            dqmDirectory = cms.string('/export/harvested/qcdSum/ahMuTauAnalyzer_woBtag')
        ##),
        ##TTplusJets = processEntryAHtoMuTau.clone(
        ##    dqmDirectory = cms.string('/export/harvested/TTplusJets/ahMuTauAnalyzer_woBtag')
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
