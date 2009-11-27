#include "TauAnalysis/FittingTools/interface/TF1WrapperBase.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <TPRegexp.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>

#include <vector>
#include <string>

TF1WrapperBase::TF1WrapperBase(const edm::ParameterSet& cfg)
  : tf1_(0),
    cfgError_(0)
{}

TF1WrapperBase::~TF1WrapperBase()
{
  delete tf1_;
}

void TF1WrapperBase::setTF1Parameter(const edm::ParameterSet& cfg)
{
  TPRegexp regexpParser_parName("par([[:digit:]]+)");

  typedef std::vector<std::string> vstring;

  vstring parNames = cfg.getParameterNamesForType<edm::ParameterSet>();
  for ( vstring::const_iterator parName = parNames.begin();
	parName != parNames.end(); ++parName ) {
    edm::ParameterSet cfgParameter = cfg.getParameter<edm::ParameterSet>(*parName);

    TString parName_tstring = parName->data();

//--- check if parameter name matches format "par0".."parN"
    if ( !regexpParser_parName.Match(parName_tstring) == 1 ) {
      edm::LogError("TF1WrapperBase::setTF1Parameter") << " Failed to parse parameter Name = " << (*parName) << " !!";
      cfgError_ = 1;
      continue;
    }

//--- extract parameter identifier;
//    (the matched integer number is referred to by the second match;
//     the first match refers to the entire parameter name)
    TObjArray* subStrings = regexpParser_parName.MatchS(parName_tstring);
    int parId = ((TObjString*)subStrings->At(1))->GetString().Atoi();
    //std::cout << "parName = " << (*parName) << ": parId = " << parId << std::endl;

    double parValue_initial = cfgParameter.getParameter<double>("initial");
    tf1_->SetParameter(parId, parValue_initial);
    
    if ( cfgParameter.exists("min") &&
	 cfgParameter.exists("max") ) {
      double parValue_min = cfgParameter.getParameter<double>("min");
      double parValue_max = cfgParameter.getParameter<double>("max");
      tf1_->SetParLimits(parId, parValue_min, parValue_max);
    }

    tf1InitialValues_[parId] = parValue_initial;
  }
}

void TF1WrapperBase::reinitializeTF1Parameter()
{
  for ( std::map<int, double>::const_iterator tf1InitialValue = tf1InitialValues_.begin();
	tf1InitialValue != tf1InitialValues_.end(); ++tf1InitialValue ) {
    tf1_->SetParameter(tf1InitialValue->first, tf1InitialValue->second);
  }
}

#include "FWCore/Framework/interface/MakerMacros.h"

EDM_REGISTER_PLUGINFACTORY(TF1WrapperPluginFactory, "TF1WrapperPluginFactory");


