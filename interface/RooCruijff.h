/*****************************************************************************
 * Package: RooRarFit                                                        *
 *    File: $Id: RooCruijff.rdl,v 1.2 2007/06/29 08:37:34 zhanglei Exp $   *
 * Authors:                                                                  *
 *    Karsten Koeneke, Massachusetts Institute of Technology, Cambridge, USA *
 *    Vouter Hulbergen                                                       *
 *                                                                           *
 * Copyright (c) 2006, Massachsetts Institute of Technology, Cambridge, USA  *
 *****************************************************************************/

#ifndef ROO_CRUIJFF
#define ROO_CRUIJFF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "TauAnalysis/FittingTools/interface/RooAbsEstimatablePdf.h"

class RooRealVar;
class RooAbsReal;

class RooCruijff : public RooAbsPdf {
public:
  RooCruijff():RooAbsPdf(){}

  RooCruijff(const char *name, const char *title,
	     RooAbsReal& _x,
	     RooAbsReal& _m0,
	     RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
	     RooAbsReal& _alphaL, RooAbsReal& _alphaR);

  RooCruijff(const RooCruijff& other, const char* name = 0);

  virtual TObject* clone(const char* newname) const {
    return new RooCruijff(*this,newname); }

  virtual RooArgSet estimateParameters(
      RooAbsData& data, double errorFactor=1.0);

  inline virtual ~RooCruijff() { }

protected:
  RooRealProxy x;
  RooRealProxy m0;
  RooRealProxy sigmaL;
  RooRealProxy sigmaR;
  RooRealProxy alphaL;
  RooRealProxy alphaR;

  Double_t evaluate() const;

private:
  ClassDef(RooCruijff,0)
};

#endif

