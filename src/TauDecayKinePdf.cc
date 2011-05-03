#include "TauAnalysis/FittingTools/interface/TauDecayKinePdf.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>

#include <TMath.h>

#include <string>
#include <iomanip>

TauDecayKinePdf::TauDecayKinePdf()
  : RooAbsPdf(),
    doAnalyticIntegration_(true),
    verbosity_(0)
{
  if ( this->verbosity_ ) std::cout << "<TauDecayKinePdf::TauDecayKinePdf(1)>:" << std::endl;
}

TauDecayKinePdf::TauDecayKinePdf(
  const char* name, const char* title, 
  RooAbsReal& x, 
  RooAbsReal& gmean, RooAbsReal& gsigma, RooAbsReal& slope, RooAbsReal& offset,  RooAbsReal& C,
  RooAbsReal& mp1, RooAbsReal& width1, RooAbsReal& mp2, RooAbsReal& width2, 
  RooAbsReal& x0, RooAbsReal& dx1)
  : RooAbsPdf(name, title), 
    x_("x", "x", this, x),
    gmean_("gmean", "gmean", this, gmean),
    gsigma_("gsigma", "gsigma", this, gsigma),
    slope_("slope", "slope", this, slope),
    offset_("offset", "offset", this, offset),
    C_("C", "C", this, C),
    mp1_("mp1", "mp1", this, mp1),
    width1_("width1", "width1", this, width1),
    mp2_("mp2", "mp2", this, mp2),
    width2_("width2", "width2", this, width2),
    x0_("x0", "x0", this, x0),
    dx1_("dx1", "dx1", this, dx1),
    doAnalyticIntegration_(true),
    verbosity_(0)
{
  if ( this->verbosity_ ) {
    std::cout << "<TauDecayKinePdf::TauDecayKinePdf(2)>:" << std::endl;
    this->print(std::cout);
  }
}

TauDecayKinePdf::TauDecayKinePdf(const TauDecayKinePdf& bluePrint, const char* newName)
  : RooAbsPdf(bluePrint, newName), 
    x_("x", this, bluePrint.x_),
    gmean_("gmean", this, bluePrint.gmean_),
    gsigma_("gsigma", this, bluePrint.gsigma_),
    slope_("slope", this, bluePrint.slope_),
    offset_("offset", this, bluePrint.offset_),
    C_("C", this, bluePrint.C_),
    mp1_("mp1", this, bluePrint.mp1_),
    width1_("width1", this, bluePrint.width1_),
    mp2_("mp2", this, bluePrint.mp2_),
    width2_("width2", this, bluePrint.width2_),
    x0_("x0", this, bluePrint.x0_),
    dx1_("dx1", this, bluePrint.dx1_),
    doAnalyticIntegration_(bluePrint.doAnalyticIntegration_),
    verbosity_(0)
{
  if ( this->verbosity_ ) {
    std::cout << "<TauDecayKinePdf::TauDecayKinePdf(3)>:" << std::endl;
    this->print(std::cout);
  }
}
  
TauDecayKinePdf::~TauDecayKinePdf()
{
//--- nothing to be done yet...
}

Double_t TauDecayKinePdf::evaluate() const
{
  if ( this->verbosity_ ) {
    std::cout << "<TauDecayKinePdf::evaluate>:" << std::endl;
    this->print(std::cout);
  }

  Double_t retVal = 0.;

  if      ( x_ <  x0_         ) retVal = evaluateGaussian(x_);
  else if ( x_ < (x0_ + dx1_) ) retVal = evaluateLandau(x_);
  else                          retVal = evaluateExponential(x_);
  
  assert(retVal >= 0.);

  return retVal;
}

Double_t TauDecayKinePdf::evaluateGaussian(Double_t x) const
{
  double gaussian = 0.;
  if ( gsigma_ > 0. ) {
    Double_t pull = (x - gmean_)/gsigma_;
    gaussian += C_*TMath::Exp(-0.5*pull*pull);
  }
  
  gaussian += (1. - C_)*(slope_*x_ + offset_);

  if ( this->verbosity_ ) std::cout << "--> returning gaussian = " << gaussian << std::endl;
  return gaussian;
}

Double_t TauDecayKinePdf::evaluateLandau(Double_t x, bool updateNormFactor) const
{
//--- normalize Landau such that 
//      exp(-((x - gmean)/g2sigma)^2) + C*exp(-((x - mean)/g4sigma)^4) = Landau(x, mp, width) 
//    at x = x0
  if ( updateNormFactor ) updateNormFactor0to1();
  double landau = ( width1_ > 0. ) ? (norm0to1_*::ROOT::Math::landau_pdf((x - mp1_)/width1_)) : 0.;
  if ( this->verbosity_ ) std::cout << "--> returning landau = " << landau << std::endl;
  return landau;
}

Double_t TauDecayKinePdf::evaluateExponential(Double_t x) const
{
//--- normalize exponential such that
//       Landau(x, mp1, width1) = Exponential(x, mp2, width2)
//    at x = x1       
  updateNormFactor0to1();
  updateNormFactor1to2();
  double exponential = ( width2_ > 0. ) ? (norm1to2_*::ROOT::Math::landau_pdf((x - mp2_)/width2_)) : 0.;
  if ( this->verbosity_ ) std::cout << "--> returning exponential = " << exponential << std::endl;
  return exponential;  
}

void TauDecayKinePdf::updateNormFactor0to1() const
{
  if ( this->verbosity_ ) std::cout << "<updateNormFactor0to1>:" << std::endl;
  double landau1_x0 = ( width1_ > 0. ) ? (::ROOT::Math::landau_pdf((x0_ - mp1_)/width1_)) : -1.;
  norm0to1_ = ( landau1_x0 > 0. ) ? (evaluateGaussian(x0_)/landau1_x0) : 0.;
  if ( this->verbosity_ ) std::cout << "--> setting norm0to1 = " << norm0to1_ << std::endl;
}

void TauDecayKinePdf::updateNormFactor1to2() const
{
  if ( this->verbosity_ ) std::cout << "<updateNormFactor1to2>:" << std::endl;
  Double_t x1 = x0_ + dx1_;
  Double_t landau2_x1 = ( width2_ > 0. ) ? (::ROOT::Math::landau_pdf((x1 - mp2_)/width2_)) : -1.;
  norm1to2_ = ( landau2_x1 > 0. ) ? (evaluateLandau(x1, false)/landau2_x1) : 0.;
  if ( this->verbosity_ ) std::cout << "--> setting norm1to2 = " << norm1to2_ << std::endl;
}

//
//-------------------------------------------------------------------------------
//

Int_t TauDecayKinePdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char*) const 
{
  if ( doAnalyticIntegration_ && matchArgs(allVars, analVars, x_) ) return 1;
  return 0 ;
}

Double_t TauDecayKinePdf::analyticalIntegral(Int_t code, const char* rangeName) const 
{
  if ( this->verbosity_ ) std::cout << "<TauDecayKinePdf::analyticalIntegral>:" << std::endl;
  
  assert(code == 1);

  Double_t retVal = 0.;

  if ( code == 1 ) {  
    updateNormFactor0to1();
    updateNormFactor1to2();

    double xMin = x_.min(rangeName);
    double xMax = x_.max(rangeName);
    
    double gaussian_integral = 0.;

    double xMin_gaussian = xMin;
    double xMax_gaussian = TMath::Min(xMax, x0_);

    if ( xMax_gaussian > xMin_gaussian ) {
      if ( gsigma_ > 0. ) {
	Double_t sqrt2_times_sigma = TMath::Sqrt(2.)*gsigma_;
	gaussian_integral += 
	  C_*TMath::Sqrt(0.5*TMath::Pi())*gsigma_
	 *(TMath::Erf((xMax_gaussian - gmean_)/sqrt2_times_sigma) - TMath::Erf((xMin_gaussian - gmean_)/sqrt2_times_sigma));
      }

      gaussian_integral += 
	(1. - C_)*(0.5*slope_*(xMax_gaussian*xMax_gaussian - xMin_gaussian*xMin_gaussian)  + offset_*(xMax_gaussian - xMin_gaussian));
    }

    if ( this->verbosity_ ) std::cout << "--> gaussian_integral = " << gaussian_integral << std::endl;
    
    double landau_integral = 0.;

    double xMin_landau = TMath::Max(xMin, x0_);
    double xMax_landau = TMath::Min(xMax, x0_ + dx1_);
    
    if ( xMax_landau > xMin_landau ) {
      landau_integral += 
	norm0to1_*width1_
       *(::ROOT::Math::landau_cdf((xMax_landau - mp1_)/width1_) - ::ROOT::Math::landau_cdf((xMin_landau - mp1_)/width1_));
    }
     
    if ( this->verbosity_ ) std::cout << "--> landau_integral = " << landau_integral << std::endl;
    
    double exponential_integral = 0.;

    double xMin_exponential = TMath::Max(xMin, x0_ + dx1_);
    double xMax_exponential = xMax;

    if ( xMax_exponential > xMin_exponential ) {
      exponential_integral += 
        norm1to2_*width2_
       *(::ROOT::Math::landau_cdf((xMax_exponential - mp2_)/width2_) - ::ROOT::Math::landau_cdf((xMin_exponential - mp2_)/width2_));
    }

    if ( this->verbosity_ ) std::cout << "--> exponential_integral = " << exponential_integral << std::endl;
    
    retVal = gaussian_integral + landau_integral + exponential_integral;
  }

  return retVal;
}
 
//
//-------------------------------------------------------------------------------
//

void TauDecayKinePdf::print(std::ostream& stream) const
{
  stream << "<TauDecayKinePdf::print>:" << std::endl;
  stream << " x: name = " << x_.absArg()->GetName() << ", value = " << x_ << std::endl;
  stream << " gmean: name = " << gmean_.absArg()->GetName() << ", value = " << gmean_ << std::endl;
  stream << " gsigma: name = " << gsigma_.absArg()->GetName() << ", value = " << gsigma_ << std::endl;
  stream << " slope: name = " << slope_.absArg()->GetName() << ", value = " << slope_ << std::endl;
  stream << " offset: name = " << offset_.absArg()->GetName() << ", value = " << offset_ << std::endl;
  stream << " C: name = " << C_.absArg()->GetName() << ", value = " << C_ << std::endl;
  stream << " mp1: name = " << mp1_.absArg()->GetName() << ", value = " << mp1_ << std::endl;
  stream << " width1: name = " << width1_.absArg()->GetName() << ", value = " << width1_ << std::endl;
  stream << " mp2: name = " << mp2_.absArg()->GetName() << ", value = " << mp2_ << std::endl;
  stream << " width2: name = " << width2_.absArg()->GetName() << ", value = " << width2_ << std::endl;
  stream << " x0: name = " << x0_.absArg()->GetName() << ", value = " << x0_ << std::endl;
  stream << " dx1: name = " << dx1_.absArg()->GetName() << ", value = " << dx1_ << std::endl;
}

// generate CInt dictionaries
ClassImp(TauDecayKinePdf)

