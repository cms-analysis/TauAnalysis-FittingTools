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
  RooAbsReal& gmean, RooAbsReal& gsigma, RooAbsReal& slope, RooAbsReal& offset,  RooAbsReal& C,
  //RooAbsReal& kappa, RooAbsReal& theta, RooAbsReal& mu, 
  RooAbsReal& mp, RooAbsReal& width, RooAbsReal& alpha,
  RooAbsReal& x0, RooAbsReal& dx1)
  : RooAbsPdf(name, title), 
    x_("x", "x", this, x),
    gmean_("gmean", "gmean", this, gmean),
    gsigma_("gsigma", "gsigma", this, gsigma),
    slope_("slope", "slope", this, slope),
    offset_("offset", "offset", this, offset),
    C_("C", "C", this, C),
    //kappa_("kappa", "kappa", this, kappa),
    //theta_("theta", "theta", this, theta),
    //mu_("mu", "mu", this, mu),
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
    gsigma_("gsigma", this, bluePrint.gsigma_),
    slope_("slope", this, bluePrint.slope_),
    offset_("offset", this, bluePrint.offset_),
    C_("C", this, bluePrint.C_),
    //kappa_("kappa", this, bluePrint.kappa_),
    //theta_("theta", this, bluePrint.theta_),
    //mu_("mu", this, bluePrint.mu_),
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
  //if      ( x_ <  x0_         ) retVal = evaluateExpGamma(x_);
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
    //gaussian = TMath::Exp(-0.5*pull*pull);
  }
  
  gaussian += (1. - C_)*(slope_*x_ + offset_);

  std::cout << "--> returning gaussian = " << gaussian << std::endl;
  return gaussian;
}
/*
Double_t TauDecayKinePdf::evaluateExpGamma(Double_t x) const
{
  double expGamma = 0.;

  if ( theta_ > 0. ) {
    double pull = (x - mu_)/theta_;
    expGamma = (1./(theta_*TMath::Gamma(kappa_)))*TMath::Exp(-TMath::Exp(pull) + kappa_*pull);
  }

  //std::cout << "--> returning expGamma = " << expGamma << std::endl;
  return expGamma;
}
 */
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
  //std::cout << "<updateNormFactor0to1>:" << std::endl;
  double landau_x0 = ( width_ > 0. ) ? (::ROOT::Math::landau_pdf((x0_ - mp_)/width_)) : -1.;
  norm0to1_ = ( landau_x0 > 0. ) ? (evaluateGaussian(x0_)/landau_x0) : 0.;
  //norm0to1_ = ( landau_x0 > 0. ) ? (evaluateExpGamma(x0_)/landau_x0) : 0.;
  //std::cout << "--> setting norm0to1 = " << norm0to1_ << std::endl;
}

void TauDecayKinePdf::updateNormFactor1to2() const
{
  //std::cout << "<updateNormFactor1to2>:" << std::endl;
  Double_t x1 = x0_ + dx1_;
  Double_t landau_x1 = ( width_ > 0. ) ? (::ROOT::Math::landau_pdf((x1 - mp_)/width_)) : -1.;
  norm1to2_ = ( landau_x1 > 0. ) ? (norm0to1_*landau_x1) : 0.;
  //std::cout << "--> setting norm1to2 = " << norm1to2_ << std::endl;
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
  //std::cout << "<TauDecayKinePdf::analyticalIntegral>:" << std::endl;
  
  assert(code == 1);

  Double_t retVal = 0.;

  if ( code == 1 ) {  
    updateNormFactor0to1();
    updateNormFactor1to2();
    
    double gaussian_integral = 0.;

    if ( g2sigma_ > 0. ) {
      Double_t sqrt2_times_sigma = TMath::Sqrt(2.)*g2sigma_;
      gaussian_integral += (1. - C_)*TMath::Sqrt(0.5*TMath::Pi())*g2sigma_
               *(TMath::Erf(gmean_/sqrt2_times_sigma) + TMath::Erf((x0_ - gmean_)/sqrt2_times_sigma));
    }

    if ( g4sigma_ > 0. ) {
      Double_t pull4 = (x0_ - gmean_)/g4sigma_;
      gaussian_integral += C_*TMath::Power(2., 1.75)*g4sigma_*(TMath::Gamma(0.25, 0.5*pull4*pull4) - TMath::Gamma(0.25, 0.));
    }

    //std::cout << "--> gaussian_integral = " << gaussian_integral << std::endl;

    double x1 = x0_ + dx1_;
    if ( x1 > x_.max(rangeName) ) x1 = x_.max(rangeName);

    double landau_integral = norm0to1_*width_
                            *(::ROOT::Math::landau_cdf((x1 - mp_)/width_) - ::ROOT::Math::landau_cdf((x0_ - mp_)/width_));
    //std::cout << "--> landau_integral = " << landau_integral << std::endl;
    
    double exponential_integral = norm1to2_*(1./alpha_)*TMath::Exp(-alpha_*x1) - TMath::Exp(-alpha_*x_.max(rangeName));
    //std::cout << "--> exponential_integral = " << exponential_integral << std::endl;
    
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
  //stream << " kappa: name = " << kappa_.absArg()->GetName() << ", value = " << kappa_ << std::endl;
  //stream << " theta: name = " << theta_.absArg()->GetName() << ", value = " << theta_ << std::endl;
  //stream << " mu: name = " << mu_.absArg()->GetName() << ", value = " << mu_ << std::endl;
  stream << " mp: name = " << mp_.absArg()->GetName() << ", value = " << mp_ << std::endl;
  stream << " width: name = " << width_.absArg()->GetName() << ", value = " << width_ << std::endl;
  stream << " alpha: name = " << alpha_.absArg()->GetName() << ", value = " << alpha_ << std::endl;
  stream << " x0: name = " << x0_.absArg()->GetName() << ", value = " << x0_ << std::endl;
  stream << " dx1: name = " << dx1_.absArg()->GetName() << ", value = " << dx1_ << std::endl;
}
