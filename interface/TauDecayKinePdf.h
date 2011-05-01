#ifndef TauAnalysis_FittingTools_TauDecayKinePdf_h
#define TauAnalysis_FittingTools_TauDecayKinePdf_h

/** \class TauDecayKinePdf
 *
 * Class to parametrize probability density for tau leptons of given energy/Pt
 * to decay into systems of visible + invisible decay products (neutrinos produced in tau decay)
 * with a certain angle/eta-phi distance between the visible and invisible momenta.
 *
 * The probability density is modelled by a piece-wise continuous function, consisting of three pieces:
 *
 *  (1) exp(-((x-gmean)/g2sigma)^2) + C*exp(-((x-mean)/g4sigma)^4) if x <= x0 with gmean > 0, gNsigma > 0, C in [0..1]
 *  (2) Landau(x, mp, width)                                       if d0 < x <= x1 with mp real, width > 0
 *  (3) exp(-alpha*x)                                              if      x >  x1 with alpha > 0
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.4 $
 *
 * $Id: TauDecayKinePdf.h,v 1.4 2011/04/30 18:55:54 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

#include <ostream>

class TauDecayKinePdf : public RooAbsPdf 
{
 public:
  // constructor(s)
  TauDecayKinePdf();   
  TauDecayKinePdf(const char*, const char*, 
		  RooAbsReal&, RooAbsReal&, RooAbsReal&, RooAbsReal&, RooAbsReal&, RooAbsReal&, 
		  RooAbsReal&, RooAbsReal&, RooAbsReal&, RooAbsReal&, 
		  RooAbsReal&, RooAbsReal&);
  TauDecayKinePdf(const TauDecayKinePdf&, const char* = "");
  
  // destructor
  virtual ~TauDecayKinePdf();

  virtual TObject* clone(const char* newName) const { return new TauDecayKinePdf(*this, newName); }

  Int_t getAnalyticalIntegral(RooArgSet&, RooArgSet&, const char* = 0) const;
  Double_t analyticalIntegral(Int_t, const char* = 0) const;
 
  void print(std::ostream&) const;

 protected:
  Double_t evaluate() const;

  Double_t evaluateGaussian(Double_t) const;
  Double_t evaluateLandau(Double_t, bool = true) const;
  Double_t evaluateExponential(Double_t) const;

  void updateNormFactor0to1() const;
  void updateNormFactor1to2() const;

  RooRealProxy x_;

  RooRealProxy gmean_;
  RooRealProxy gsigma_;
  RooRealProxy slope_;
  RooRealProxy offset_;
  RooRealProxy C_;
  RooRealProxy mp1_;
  RooRealProxy width1_;
  RooRealProxy mp2_;
  RooRealProxy width2_;
  RooRealProxy x0_;
  RooRealProxy dx1_;

  mutable Double_t norm0to1_;
  mutable Double_t norm1to2_;

  int verbosity_;
};

#endif  
