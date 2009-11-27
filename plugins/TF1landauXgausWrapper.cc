#include "TauAnalysis/FittingTools/plugins/TF1landauXgausWrapper.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <TMath.h>

#include <string>

Double_t langaufun(Double_t* x, Double_t* par) 
{
  //-----------------------------------------------------------------------------
  // Convoluted Landau and Gaussian Fitting Function
  //         (using ROOT's Landau and Gauss functions)
  //
  //  Based on a Fortran code by R.Fruehwirth (fruhwirth@hephy.oeaw.ac.at)
  //  Adapted for C++/ROOT by H.Pernegger (Heinz.Pernegger@cern.ch) and
  //   Markus Friedl (Markus.Friedl@cern.ch)
  //
  // Downloaded from http://root.cern.ch/root/roottalk/roottalk02/att-3361/01-langaus.C
  //
  //-----------------------------------------------------------------------------

  // function parameter
  // ------------------
  //  o par[0]: Width (scale) parameter of Landau density
  //  o par[1]: Most Probable (MP, location) parameter of Landau density
  //  o par[2]: Total area (integral -inf to inf, normalization constant)
  //  o par[3]: Width (sigma) of convoluted Gaussian function
  //
  // In the Landau distribution (represented by the CERNLIB approximation), 
  // the maximum is located at x = -0.22278298 with the location parameter = 0.
  // This shift is corrected within this function, so that the actual
  // maximum is identical to the MP parameter.

  // Numeric constants
  const Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
  const Double_t mpshift  = -0.22278298;     // Landau maximum location

  // Control constants
  const Double_t np = 100.0; // number of convolution steps
  const Double_t sc =   5.0; // convolution extends to +-sc Gaussian sigmas

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow, xupp;
  Double_t step;
  Double_t i;

  // MP shift correction
  mpc = par[1] - mpshift * par[0]; 

  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];

  step = (xupp - xlow) / np;

  // Approximate convolution integral of Landau and Gaussian by sum
  for( i = 1.0; i <= np/2; i++ ) {
    xx = xlow + (i - 0.5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
    
    xx = xupp - (i - 0.5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
  }

  return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1landauXgausWrapper::TF1landauXgausWrapper(const edm::ParameterSet& cfg)
  : TF1WrapperBase(cfg) 
{
  std::string pluginName = ( cfg.exists("pluginName") ) ? cfg.getParameter<std::string>("pluginName") : "TF1landauXgausWrapper::tf1_";
  
  double xMin = cfg.getParameter<double>("xMin");
  double xMax = cfg.getParameter<double>("xMax");
  
  tf1_ = new TF1(pluginName.data(), langaufun, xMin, xMax, 4);
  
  tf1_->SetParNames("Width", "MP", "Area", "GSigma");

  setTF1Parameter(cfg.getParameter<edm::ParameterSet>("parameter"));
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(TF1WrapperPluginFactory, TF1landauXgausWrapper, "TF1landauXgausWrapper");
