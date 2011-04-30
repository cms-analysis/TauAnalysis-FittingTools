#include "TauAnalysis/FittingTools/interface/TauDecayKinePdf.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>

#include <TMath.h>

#include <string>
#include <iomanip>

TauDecayKinePdf::TauDecayKinePdf()
  : verbosity_(0)
{
  if ( this->verbosity_ ) std::cout << "<TauDecayKinePdf::TauDecayKinePdf(1)>:" << std::endl;
}

TauDecayKinePdf::TauDecayKinePdf(
  const char* name, const char* title, 
  RooAbsReal& x, 
  RooAbsReal& gmean, RooAbsReal& gsigma, RooAbsReal& slope, RooAbsReal& offset,  RooAbsReal& C,
  RooAbsReal& mp, RooAbsReal& width, RooAbsReal& alpha,
  RooAbsReal& x0, RooAbsReal& dx1)
  : RooAbsPdf(name, title), 
    x_("x", "x", this, x),
    gmean_("gmean", "gmean", this, gmean),
    gsigma_("gsigma", "gsigma", this, gsigma),
    slope_("slope", "slope", this, slope),
    offset_("offset", "offset", this, offset),
    C_("C", "C", this, C),
    mp_("mp", "mp", this, mp),
    width_("width", "width", this, width),
    alpha_("alpha", "alpha", this, alpha),
    x0_("x0", "x0", this, x0),
    dx1_("dx1", "dx1", this, dx1),
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
    mp_("mp", this, bluePrint.mp_),
    width_("width", this, bluePrint.width_),
    alpha_("alpha", this, bluePrint.alpha_),
    x0_("x0", this, bluePrint.x0_),
    dx1_("dx1", this, bluePrint.dx1_),
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

Double_t TauDecayKinePdf::evaluateLandau(Double_t x) const
{
//--- normalize Landau such that 
//      exp(-((x - gmean)/g2sigma)^2) + C*exp(-((x - mean)/g4sigma)^4) = Landau(x, mp, width) 
//    at x = x0
  updateNormFactor0to1();
  double landau = ( width_ > 0. ) ? (norm0to1_*::ROOT::Math::landau_pdf((x - mp_)/width_)) : 0.;
  if ( this->verbosity_ ) std::cout << "--> returning landau = " << landau << std::endl;
  return landau;
}

Double_t TauDecayKinePdf::evaluateExponential(Double_t x) const
{
//--- normalize exponential such that
//       Landau(x, mp, width) = exp(-alpha*x) 
//    at x = x1       
  updateNormFactor0to1();
  updateNormFactor1to2();
  Double_t x1 = x0_ + dx1_;
  double exponential = norm1to2_*TMath::Exp(-alpha_*(x - x1));
  if ( this->verbosity_ ) std::cout << "--> returning exponential = " << exponential << std::endl;
  return exponential;  
}

void TauDecayKinePdf::updateNormFactor0to1() const
{
  if ( this->verbosity_ ) std::cout << "<updateNormFactor0to1>:" << std::endl;
  double landau_x0 = ( width_ > 0. ) ? (::ROOT::Math::landau_pdf((x0_ - mp_)/width_)) : -1.;
  norm0to1_ = ( landau_x0 > 0. ) ? (evaluateGaussian(x0_)/landau_x0) : 0.;
  if ( this->verbosity_ ) std::cout << "--> setting norm0to1 = " << norm0to1_ << std::endl;
}

void TauDecayKinePdf::updateNormFactor1to2() const
{
  if ( this->verbosity_ ) std::cout << "<updateNormFactor1to2>:" << std::endl;
  Double_t x1 = x0_ + dx1_;
  Double_t landau_x1 = ( width_ > 0. ) ? (::ROOT::Math::landau_pdf((x1 - mp_)/width_)) : -1.;
  norm1to2_ = ( landau_x1 > 0. ) ? (norm0to1_*landau_x1) : 0.;
  if ( this->verbosity_ ) std::cout << "--> setting norm1to2 = " << norm1to2_ << std::endl;
}

//
//-------------------------------------------------------------------------------
//
/*
Int_t TauDecayKinePdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char*) const 
{
  if ( matchArgs(allVars, analVars, x_) ) return 1 ;
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
    
    double gaussian_integral = 0.;

    if ( gsigma_ > 0. ) {
      Double_t sqrt2_times_sigma = TMath::Sqrt(2.)*gsigma_;
      gaussian_integral += C_*TMath::Sqrt(0.5*TMath::Pi())*gsigma_
               *(TMath::Erf(gmean_/sqrt2_times_sigma) + TMath::Erf((x0_ - gmean_)/sqrt2_times_sigma));
    }

    gaussian_integral += (1. - C_)*(0.5*slope_*x0_*x0_ + offset_*x0_);

    if ( this->verbosity_ ) std::cout << "--> gaussian_integral = " << gaussian_integral << std::endl;

    double x1 = x0_ + dx1_;
    if ( x1 > x_.max(rangeName) ) x1 = x_.max(rangeName);

    double landau_integral = norm0to1_*width_
                            *(::ROOT::Math::landau_cdf((x1 - mp_)/width_) - ::ROOT::Math::landau_cdf((x0_ - mp_)/width_));
    if ( this->verbosity_ ) std::cout << "--> landau_integral = " << landau_integral << std::endl;
    
    double exponential_integral = norm1to2_*(1./alpha_)*TMath::Exp(-alpha_*x1) - TMath::Exp(-alpha_*x_.max(rangeName));
    if ( this->verbosity_ ) std::cout << "--> exponential_integral = " << exponential_integral << std::endl;
    
    retVal = gaussian_integral + landau_integral + exponential_integral;
  }

  return retVal;
}
 */ 
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
  stream << " mp: name = " << mp_.absArg()->GetName() << ", value = " << mp_ << std::endl;
  stream << " width: name = " << width_.absArg()->GetName() << ", value = " << width_ << std::endl;
  stream << " alpha: name = " << alpha_.absArg()->GetName() << ", value = " << alpha_ << std::endl;
  stream << " x0: name = " << x0_.absArg()->GetName() << ", value = " << x0_ << std::endl;
  stream << " dx1: name = " << dx1_.absArg()->GetName() << ", value = " << dx1_ << std::endl;
}
