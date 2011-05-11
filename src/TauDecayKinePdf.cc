#include "TauAnalysis/FittingTools/interface/TauDecayKinePdf.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//#include "TauAnalysis/FittingTools/interface/owen_T.h"
#include "TauAnalysis/FittingTools/src/owens.H"

#include "RooCmdArg.h"

#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>

#include <TMath.h>
#include <TH1.h>

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
  RooAbsReal& gmean, RooAbsReal& gsigma, RooAbsReal& alpha, RooAbsReal& slope, RooAbsReal& offset,  RooAbsReal& C,
  RooAbsReal& mp1, RooAbsReal& width1, RooAbsReal& mp2, RooAbsReal& width2,
  RooAbsReal& x0, RooAbsReal& dx1)
  : RooAbsPdf(name, title),
    x_("x", "x", this, x),
    gmean_("gmean", "gmean", this, gmean),
    gsigma_("gsigma", "gsigma", this, gsigma),
    alpha_("alpha", "alpha", this, alpha),
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
    alpha_("alpha", this, bluePrint.alpha_),
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

  x0Value_  = TMath::Max(x0_,  0.);
  dx1Value_ = TMath::Max(dx1_, 0.);

  //-----------------------------------------------------------------------------
  // CV: "hack" to resolve ambiguity between parameters
  //     of first and second Landau distribution in case x0 + dx1 > 11.5
  //
  //    (second Landau distribution models tails at large x,
  //     where event statistics is low and width not well constrained)
  //
  //    --> take width from fit of first Landau distribution instead
  //
  if ( dx1_ < 1.e-6 ) {
    mp2Value_ = mp1_;
    width2Value_ = width1_;
  } else {
    mp2Value_ = mp2_;
    width2Value_ = width2_;
  }
  //-----------------------------------------------------------------------------

  Double_t retVal = 0.;

  if      ( x_ <  x0Value_              ) retVal = evaluateGaussian(x_);
  else if ( x_ < (x0Value_ + dx1Value_) ) retVal = evaluateLandau1(x_);
  else                                    retVal = evaluateLandau2(x_);

  const double epsilon = 1.e-12;
  if ( retVal < epsilon || TMath::IsNaN(retVal) ) retVal = epsilon;

  return retVal;
}

Double_t TauDecayKinePdf::evaluateGaussian(Double_t x) const
{
  double gaussian = 0.;
  if ( gsigma_ > 0. ) {
    Double_t pull = (x - gmean_)/gsigma_;
    gaussian += C_*(1./(TMath::Sqrt(2.*TMath::Pi())*gsigma_))*TMath::Exp(-0.5*pull*pull)*(1. + TMath::Erf(alpha_*pull));
  }

  gaussian += (1. - C_)*(slope_*x + offset_);

  if ( this->verbosity_ ) std::cout << "--> returning Gaussian = " << gaussian << std::endl;
  return gaussian;
}

Double_t TauDecayKinePdf::evaluateLandau1(Double_t x, bool updateNormFactor) const
{
//--- normalize Landau such that
//      exp(-((x - gmean)/g2sigma)^2) + C*exp(-((x - mean)/g4sigma)^4) = Landau1(x, mp, width)
//    at x = x0
  if ( updateNormFactor ) updateNormFactor0to1();
  double landau1 = ( width1_ > 0. ) ? (norm0to1_*::ROOT::Math::landau_pdf((x - mp1_)/width1_)) : 0.;
  if ( this->verbosity_ ) std::cout << "--> returning Landau1 = " << landau1 << std::endl;
  return landau1;
}

Double_t TauDecayKinePdf::evaluateLandau2(Double_t x) const
{
//--- normalize exponential such that
//       Landau1(x, mp1, width1) = Landau2(x, mp2, width2)
//    at x = x1
  updateNormFactor0to1();
  updateNormFactor1to2();
  double landau2 = ( width2Value_ > 0. ) ? (norm1to2_*::ROOT::Math::landau_pdf((x - mp2Value_)/width2Value_)) : 0.;
  if ( this->verbosity_ ) std::cout << "--> returning Landau2 = " << landau2 << std::endl;
  return landau2;
}

void TauDecayKinePdf::updateNormFactor0to1() const
{
  if ( this->verbosity_ ) std::cout << "<updateNormFactor0to1>:" << std::endl;
  double landau1_x0 = ( width1_ > 0. ) ? (::ROOT::Math::landau_pdf((x0Value_ - mp1_)/width1_)) : -1.;
  norm0to1_ = ( landau1_x0 > 0. ) ? (evaluateGaussian(x0Value_)/landau1_x0) : 0.;
  if ( this->verbosity_ ) std::cout << "--> setting norm0to1 = " << norm0to1_ << std::endl;
}

void TauDecayKinePdf::updateNormFactor1to2() const
{
  if ( this->verbosity_ ) std::cout << "<updateNormFactor1to2>:" << std::endl;
  Double_t x1 = x0Value_ + dx1Value_;
  Double_t landau2_x1 = ( width2Value_ > 0. ) ? (::ROOT::Math::landau_pdf((x1 - mp2Value_)/width2Value_)) : -1.;
  norm1to2_ = ( landau2_x1 > 0. ) ? (evaluateLandau1(x1, false)/landau2_x1) : 0.;
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
    double xMax_gaussian = TMath::Min(xMax, x0Value_);

    if ( xMax_gaussian > xMin_gaussian ) {
      if ( gsigma_ > 0. ) {
	Double_t sqrt2_times_sigma = TMath::Sqrt(2.)*gsigma_;
	double cdfHi = 0.5*(1. + TMath::Erf((xMax_gaussian - gmean_)/sqrt2_times_sigma));
	if ( alpha_ != 0. ) cdfHi -= 2.*t((xMax_gaussian - gmean_)/gsigma_, alpha_);
	double cdfLo = 0.5*(1. + TMath::Erf((xMin_gaussian - gmean_)/sqrt2_times_sigma));
	if ( alpha_ != 0. ) cdfLo -= 2.*t((xMin_gaussian - gmean_)/gsigma_, alpha_);
	gaussian_integral += C_*(cdfHi - cdfLo);
      }

      gaussian_integral +=
	(1. - C_)*(0.5*slope_*(xMax_gaussian*xMax_gaussian - xMin_gaussian*xMin_gaussian)  + offset_*(xMax_gaussian - xMin_gaussian));
    }

    if ( this->verbosity_ ) std::cout << "--> gaussian_integral = " << gaussian_integral << std::endl;

    double landau1_integral = 0.;

    double xMin_landau1 = TMath::Max(xMin, x0Value_);
    double xMax_landau1 = TMath::Min(xMax, x0Value_ + dx1Value_);

    if ( xMax_landau1 > xMin_landau1 ) {
      if ( width1_ > 0. ) {
	landau1_integral +=
	  norm0to1_*width1_
	 *(::ROOT::Math::landau_cdf((xMax_landau1 - mp1_)/width1_) - ::ROOT::Math::landau_cdf((xMin_landau1 - mp1_)/width1_));
      }
    }

    if ( this->verbosity_ ) std::cout << "--> Landau1_integral = " << landau1_integral << std::endl;

    double landau2_integral = 0.;

    double xMin_landau2 = TMath::Max(xMin, x0Value_ + dx1Value_);
    double xMax_landau2 = xMax;

    if ( xMax_landau2 > xMin_landau2 ) {
      if ( width2Value_ > 0. ) {
        landau2_integral +=
          norm1to2_*width2Value_
         *(::ROOT::Math::landau_cdf((xMax_landau2 - mp2Value_)/width2Value_)
          - ::ROOT::Math::landau_cdf((xMin_landau2 - mp2Value_)/width2Value_));
      }
    }

    if ( this->verbosity_ ) std::cout << "--> Landau2_integral = " << landau2_integral << std::endl;

    retVal = gaussian_integral + landau1_integral + landau2_integral;
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
  stream << " alpha: name = " << alpha_.absArg()->GetName() << ", value = " << alpha_ << std::endl;
  stream << " slope: name = " << slope_.absArg()->GetName() << ", value = " << slope_ << std::endl;
  stream << " offset: name = " << offset_.absArg()->GetName() << ", value = " << offset_ << std::endl;
  stream << " C: name = " << C_.absArg()->GetName() << ", value = " << C_ << std::endl;
  stream << " mp1: name = " << mp1_.absArg()->GetName() << ", value = " << mp1_ << std::endl;
  stream << " width1: name = " << width1_.absArg()->GetName() << ", value = " << width1_ << std::endl;
  stream << " mp2: name = " << mp2_.absArg()->GetName() << ", value = " << mp2Value_ << std::endl;
  stream << " width2: name = " << width2_.absArg()->GetName() << ", value = " << width2Value_ << std::endl;
  stream << " x0: name = " << x0_.absArg()->GetName() << ", value = " << x0Value_ << std::endl;
  stream << " dx1: name = " << dx1_.absArg()->GetName() << ", value = " << dx1Value_ << std::endl;
}

/*
 * Methods for estimating the parameters
 */


namespace {

double compHistogramRMSinRange(TH1* histogram, double xRef, int firstBin, int lastBin)
{
  double rms = 0.;

  double sumBinContents = 0.;

  for ( int iBin = firstBin; iBin <= lastBin; ++iBin ) {
    double binContent = histogram->GetBinContent(iBin);
    double x = histogram->GetBinCenter(iBin);

    double dx = x - xRef;
    rms += binContent*dx*dx;

    sumBinContents += binContent;
  }

  if ( sumBinContents != 0. ) rms = TMath::Sqrt(rms/sumBinContents);

  return rms;
}

double compHistogramRMSltMax(TH1* histogram)
{
  int binMax = histogram->GetMaximumBin();

  double xMax = histogram->GetBinCenter(binMax);

  return compHistogramRMSinRange(histogram, xMax, 1, binMax);
}

double compHistogramRMSgtMax(TH1* histogram)
{
  int binMax = histogram->GetMaximumBin();

  double xMax = histogram->GetBinCenter(binMax);

  int numBins = histogram->GetNbinsX();

  return compHistogramRMSinRange(histogram, xMax, binMax, numBins);
}

double compFallingEdgePos(TH1* histogram, double& errEdgePosRight, double& errEdgePosLeft)
{
  std::cout << "<compFallingEdgePos>:" << std::endl;

  int numBins = histogram->GetNbinsX();

  std::vector<double> window5derrivatives(numBins);

  for ( int iBin = 3; iBin <= (numBins - 2); ++iBin ) {
    double diffBinContent = histogram->GetBinContent(iBin + 2) - histogram->GetBinContent(iBin - 2);
    double diffBinCenter  = histogram->GetBinCenter(iBin + 2) - histogram->GetBinCenter(iBin - 2);

    window5derrivatives[iBin] = diffBinContent/diffBinCenter;
  }

  double minWindow5derrivative = 0;
  int binFallingEdge = -1.;
  double fallindEdgePos = -1.;
  for ( int iBin = 3; iBin <= (numBins - 2); ++iBin ) {
    if ( window5derrivatives[iBin] < minWindow5derrivative ) {
      binFallingEdge = iBin;
      fallindEdgePos = histogram->GetBinCenter(iBin);
      minWindow5derrivative = window5derrivatives[iBin];
    }
  }

  //binFallingEdge -= 2; // CV: "phenomenological" correction...

  if ( binFallingEdge > 0 ) {
    errEdgePosRight = TMath::Abs(histogram->GetBinContent(binFallingEdge)/minWindow5derrivative);
    std::cout << " errEdgePosRight = " << errEdgePosRight << std::endl;

    double histogramMax_y = histogram->GetBinContent(histogram->GetMaximumBin());
    errEdgePosLeft  = TMath::Abs((histogramMax_y - histogram->GetBinContent(binFallingEdge))/minWindow5derrivative);
    std::cout << " errEdgePosLeft = " << errEdgePosLeft << std::endl;
  }

  return fallindEdgePos;
}
}

RooArgSet TauDecayKinePdf::estimateParameters(RooAbsData& data, double errorFactor) {
  // Build a histogram of the observable
  std::auto_ptr<TH1> histo(data.createHistogram(x_.arg().GetName()));
  TH1* histogram = histo.get();
  assert(histogram);

  double histogramMean = histogram->GetMean();
  //if ( histogramMean > 0.15 ) histogramMean = 0.15;
  std::cout << " histogramMean = " << histogramMean << std::endl;
  double histogramMax_x = histogram->GetBinCenter(histogram->GetMaximumBin());
  //double histogramMax_y = histogram->GetBinContent(histogram->GetMaximumBin());
  //double risingEdge_slope = histogramMax_y/histogramMax_x;
  //if ( histogramMax_x > 0.15 ) histogramMax = 0.15;
  std::cout << " histogramMax_x = " << histogramMax_x << std::endl;

  // Find first bin that has a value 0.1*max above zero
  double startX = 0;
  double valueAtMax = histogram->GetBinContent(histogram->GetMaximumBin());
  for (int iBin = 1; iBin < histogram->GetNbinsX(); ++iBin) {
    if (histogram->GetBinContent(iBin) > 0.1*valueAtMax) {
      startX = histogram->GetBinCenter(iBin);
      break;
    }
  }

  std::cout << " firstNonZero_x = " << startX << std::endl;

  double errEdgePosRight = 0.;
  double errEdgePosLeft  = 0.;
  double fallingEdge_position = compFallingEdgePos(histogram, errEdgePosRight, errEdgePosLeft);
  std::cout << " fallingEdge_position = " << fallingEdge_position << std::endl;
  double histogramRMS = histogram->GetRMS();
  //if ( histogramRMS  > 0.15 ) histogramRMS  = 0.15;
  std::cout << " histogramRMS = " << histogramRMS << std::endl;
  double histogramRMSltMax = compHistogramRMSltMax(histogram);
  std::cout << " histogramRMSltMax = " << histogramRMSltMax << std::endl;
  double histogramRMSgtMax = compHistogramRMSgtMax(histogram);
  std::cout << " histogramRMSgtMax = " << histogramRMSgtMax << std::endl;
  double x0Min = 0.9*fallingEdge_position - errEdgePosLeft;
  double x0Max = 1.1*fallingEdge_position + errEdgePosRight;

  // Convert everything to RooRealVars
  RooRealVar* _x = convertProxy(x_);
  RooRealVar* _gmean = convertProxy(gmean_);
  RooRealVar* _gsigma = convertProxy(gsigma_);
  RooRealVar* _alpha = convertProxy(alpha_);
  RooRealVar* _slope = convertProxy(slope_);
  RooRealVar* _offset = convertProxy(offset_);
  RooRealVar* _C = convertProxy(C_);
  RooRealVar* _mp1 = convertProxy(mp1_);
  RooRealVar* _width1 = convertProxy(width1_);
  RooRealVar* _mp2 = convertProxy(mp2_);
  RooRealVar* _width2 = convertProxy(width2_);
  RooRealVar* _x0 = convertProxy(x0_);
  RooRealVar* _dx1 = convertProxy(dx1_);

  *_gmean = histogramMax_x;
  //gmean_.setMin(0.5*histogramMax_x);
  //gmean_.setMax(1.5*histogramMax_x);

  *_gsigma = histogramRMS;
  //gsigma_.setMin(0.25*TMath::Min(histogramRMSltMax, histogramRMS));
  //gsigma_.setMax(2.0*TMath::Min(histogramRMSltMax, histogramRMS));


  double slopeMin = 1.e-2/(histogramMax_x - startX);
  double slopeMax = 1.e+2/(histogramMax_x - startX);
  *_slope = 1./(histogramMax_x - startX);
  _slope->setMin(slopeMin);
  _slope->setMax(slopeMax);

  // Offset + slope*startX = 0
  double offsetVal = -1.0*_slope->getVal()*startX;
  *_offset = offsetVal;

  //*_slope = 0.;
  //*_offset = 0.;
  //_slope->setConstant(true);
  //_offset->setConstant(true);

  *_C = 0.25;
  //*_C = 0.0;

  *_mp1 = histogramMax_x;
  //mp1_.setMin(0.8*histogramMax_x);
  //mp1_.setMax(fallingEdge_position);

  *_width1 = 0.04*histogramRMSgtMax;
  //width1_.setMin(1.0e-4);
  //width1_.setMax(2.0*histogramRMSgtMax);

  *_mp2 = histogramMax_x;
  //mp2_.setMin(0.8*histogramMax_x);
  //mp2_.setMax(fallingEdge_position);
  //mp2_.setConstant(true);

  *_width2 = 0.04*histogramRMSgtMax;
  //width2_.setMin(1.0e-4);
  //width2_.setMax(2.0*histogramRMSgtMax);
  //width2_.setConstant(true);

  *_x0 = fallingEdge_position;

  *_dx1 = 20;
  //dx1_.setConstant(true);
  //
  std::cout << "Preforming gaussian prefit" << std::endl;
  histogram->GetXaxis()->SetRangeUser(0, fallingEdge_position);
  _x->setRange("gaussian", 0, fallingEdge_position);

  std::cout << "Estimating alpha" << std::endl;
  *_alpha = 0.0;
  double sampleSkew = histogram->GetSkewness();
  std::cout << " - initial sample skew is " << sampleSkew << std::endl;

  if (TMath::Abs(sampleSkew) > 0.80) {
    sampleSkew = 0.80*sampleSkew/TMath::Abs(sampleSkew);
  }
  double sampleSigma = histogram->GetRMS();
  double sampleMean = histogram->GetMean();

  std::cout << " - corrected sample skew is " << sampleSkew << std::endl;
  std::cout << " - sample mean is " << sampleMean << std::endl;
  std::cout << " - sample sigam is " << sampleSigma << std::endl;

  double thirdMomentToTwoThirds = TMath::Power(TMath::Abs(sampleSkew), 2.0/3.0);

  double absDelta = TMath::Sqrt(
      TMath::Pi()*0.5*(thirdMomentToTwoThirds)/(
        thirdMomentToTwoThirds + TMath::Power((4 - TMath::Pi())/2.0, 2.0/3.0)));

  double delta = absDelta*sampleSkew/TMath::Abs(sampleSkew);

  double estSkew = delta/TMath::Sqrt(1- delta*delta);
  double estScale = sampleSigma/TMath::Sqrt(1 - 2*delta*delta/TMath::Pi());
  double estLocation = sampleMean - estScale*delta*TMath::Sqrt(2/TMath::Pi());

  *_gmean = estLocation;
  *_gsigma = estScale;
  *_alpha = estSkew;

  _mp1->setConstant(true);
  _mp1->setMax(_x0->getVal());
  _width1->setConstant(true);
  _mp2->setConstant(true);
  _width2->setConstant(true);
  _x0->setConstant(true);
  _dx1->setConstant(true);
  *_C = 1;
  _C->setConstant(true);
  _slope->setConstant(true);
  _offset->setConstant(true);

  _alpha->setConstant(true);

  std::cout << " making copy of self " << std::endl;
  TauDecayKinePdf copy(*this, "tmp");

  copy.fitTo(data, RooFit::Range("gaussian"));

  std::cout << "Performing landau prefit" << std::endl;
  _x->setRange("landau1", 0.5*(histogramMax_x+ fallingEdge_position),
      fallingEdge_position + dx1_);

  _gmean->setConstant(true);
  _gsigma->setConstant(true);
  _alpha->setConstant(true);
  _slope->setConstant(true);
  _offset->setConstant(true);
  _C->setConstant(true);

  _mp1->setConstant(false);
  _width1->setConstant(false);

  copy.fitTo(data, RooFit::Range("landau1"));

  _gmean->setConstant(false);
  _gsigma->setConstant(false);
  //_alpha->setConstant(false);
  //_slope->setConstant(false);
  //_offset->setConstant(false);
  //_C->setConstant(false);
  _x0->setConstant(false);
  _alpha->setMax(0);
  _alpha->setMin(-10);

  RooArgSet output;

  output.add(gmean_.arg());
  output.add(gsigma_.arg());
  output.add(alpha_.arg());
  output.add(slope_.arg());
  output.add(offset_.arg());
  output.add(C_.arg());
  output.add(mp1_.arg());
  output.add(width1_.arg());
  output.add(mp2_.arg());
  output.add(width2_.arg());
  output.add(x0_.arg());
  output.add(dx1_.arg());

  return output;
}

// generate CInt dictionaries
ClassImp(TauDecayKinePdf)

