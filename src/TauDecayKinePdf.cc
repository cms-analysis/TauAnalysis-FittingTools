#include "TauAnalysis/FittingTools/interface/TauDecayKinePdf.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>

#include <TMath.h>

#include <string>
#include <iomanip>

TauDecayKinePdf::TauDecayKinePdf()
{
  //std::cout << "<TauDecayKinePdf::TauDecayKinePdf(1)>:" << std::endl;
}

TauDecayKinePdf::TauDecayKinePdf(
  const char* name, const char* title, 
  RooAbsReal& x, 
  RooAbsReal& gmean, RooAbsReal& g2sigma, RooAbsReal& g4sigma, RooAbsReal& C, RooAbsReal& mp, RooAbsReal& width, RooAbsReal& alpha,
  RooAbsReal& x0, RooAbsReal& dx1)
  : RooAbsPdf(name, title), 
    x_("x", "x", this, x),
    gmean_("gmean", "gmean", this, gmean),
    g2sigma_("g2sigma", "g2sigma", this, g2sigma),
    g4sigma_("g4sigma", "g4sigma", this, g4sigma),
    C_("C", "C", this, C),
    mp_("mp", "mp", this, mp),
    width_("width", "width", this, width),
    alpha_("alpha", "alpha", this, alpha),
    x0_("x0", "x0", this, x0),
    dx1_("dx1", "dx1", this, dx1)
{
  //std::cout << "<TauDecayKinePdf::TauDecayKinePdf(2)>:" << std::endl;
  //this->print(std::cout);
}

TauDecayKinePdf::TauDecayKinePdf(const TauDecayKinePdf& bluePrint, const char* newName)
  : RooAbsPdf(bluePrint, newName), 
    x_("x", this, bluePrint.x_),
    gmean_("gmean", this, bluePrint.gmean_),
    g2sigma_("g2sigma", this, bluePrint.g2sigma_),
    g4sigma_("g4sigma", this, bluePrint.g4sigma_),
    C_("C", this, bluePrint.C_),
    mp_("mp", this, bluePrint.mp_),
    width_("width", this, bluePrint.width_),
    alpha_("alpha", this, bluePrint.alpha_),
    x0_("x0", this, bluePrint.x0_),
    dx1_("dx1", this, bluePrint.dx1_)
{
  //std::cout << "<TauDecayKinePdf::TauDecayKinePdf(3)>:" << std::endl;
  //this->print(std::cout);
}
  
TauDecayKinePdf::~TauDecayKinePdf()
{
//--- nothing to be done yet...
}

Double_t TauDecayKinePdf::evaluate() const
{
  std::cout << "<TauDecayKinePdf::evaluate>:" << std::endl;
  this->print(std::cout);

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

  if ( g2sigma_ > 0. ) {
    Double_t pull2 = (x - gmean_)/g2sigma_;
    Double_t g2 = (1. - C_)*TMath::Exp(-0.5*pull2*pull2);
    std::cout << "--> adding g2 = " << g2 << std::endl;
    gaussian += g2;
  }
  
  if ( g4sigma_ > 0. ) {
    Double_t pull4 = (x - gmean_)/g4sigma_;
    Double_t g4 = C_*TMath::Exp(-0.5*pull4*pull4*pull4*pull4);
    std::cout << "--> adding g4 = " << g4 << std::endl;
    gaussian += g4;
  }
  
  std::cout << "--> returning g2 + g4 = " << gaussian << std::endl;
  return gaussian;
}

Double_t TauDecayKinePdf::evaluateLandau(Double_t x) const
{
//--- normalize Landau such that 
//      exp(-((x - gmean)/g2sigma)^2) + C*exp(-((x - mean)/g4sigma)^4) = Landau(x, mp, width) 
//    at x = x0
  updateNormFactor0to1();
  double landau = ( width_ > 0. ) ? (norm0to1_*::ROOT::Math::landau_pdf((x - mp_)/width_)) : 0.;
  std::cout << "--> returning landau = " << landau << std::endl;
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
  std::cout << "--> returning exponential = " << exponential << std::endl;
  return exponential;  
}

void TauDecayKinePdf::updateNormFactor0to1() const
{
  std::cout << "<updateNormFactor0to1>:" << std::endl;
  double landau_x0 = ( width_ > 0. ) ? (::ROOT::Math::landau_pdf((x0_ - mp_)/width_)) : -1.;
  norm0to1_ = ( landau_x0 > 0. ) ? (evaluateGaussian(x0_)/landau_x0) : 0.;
  std::cout << "--> setting norm0to1 = " << norm0to1_ << std::endl;
}

void TauDecayKinePdf::updateNormFactor1to2() const
{
  std::cout << "<updateNormFactor1to2>:" << std::endl;
  Double_t x1 = x0_ + dx1_;
  Double_t landau_x1 = ( width_ > 0. ) ? (::ROOT::Math::landau_pdf((x1 - mp_)/width_)) : -1.;
  norm1to2_ = ( landau_x1 > 0. ) ? (norm0to1_*landau_x1) : 0.;
  std::cout << "--> setting norm1to2 = " << norm1to2_ << std::endl;
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

Double_t TauDecayKinePdf::analyticalIntegral(Int_t code, const char*) const 
{
  assert(code == 1) ;

  Double_t retVal = 0.;

  if ( code == 1 ) {  
    updateNormFactor0to1();
    updateNormFactor1to2();
    
    if ( g2sigma_ > 0. ) {
      Double_t sqrt2_times_sigma = TMath::Sqrt(2.)*g2sigma_;
      retVal += TMath::Sqrt(0.5*TMath::Pi())*g2sigma_
               *(TMath::Erf(gmean_/sqrt2_times_sigma) + TMath::Erf((x0_ - gmean_)/sqrt2_times_sigma));
    }

    if ( g4sigma_ > 0. ) {
      Double_t pull4 = (x0_ - gmean_)/g4sigma_;
      retVal += TMath::Power(2., 1.75)*g4sigma_*(TMath::Gamma(0.25, 0.5*pull4*pull4) - TMath::Gamma(0.25, 0.));
    }

    double x1 = x0_ + dx1_;

    retVal += norm0to1_*(::ROOT::Math::landau_cdf(x1, width_, mp_) - ::ROOT::Math::landau_cdf(x0_, width_, mp_));

    retVal += norm1to2_*(1./alpha_)*TMath::Exp(-alpha_*x1);
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
  stream << " g2sigma: name = " << g2sigma_.absArg()->GetName() << ", value = " << g2sigma_ << std::endl;
  stream << " g4sigma: name = " << g4sigma_.absArg()->GetName() << ", value = " << g4sigma_ << std::endl;
  stream << " C: name = " << C_.absArg()->GetName() << ", value = " << C_ << std::endl;
  stream << " mp: name = " << mp_.absArg()->GetName() << ", value = " << mp_ << std::endl;
  stream << " width: name = " << width_.absArg()->GetName() << ", value = " << width_ << std::endl;
  stream << " alpha: name = " << alpha_.absArg()->GetName() << ", value = " << alpha_ << std::endl;
  stream << " x0: name = " << x0_.absArg()->GetName() << ", value = " << x0_ << std::endl;
  stream << " dx1: name = " << dx1_.absArg()->GetName() << ", value = " << dx1_ << std::endl;
}
