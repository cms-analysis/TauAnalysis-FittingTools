#ifndef TauAnalysis_FittingTools_RooSkewNormal_h
#define TauAnalysis_FittingTools_RooSkewNormal_h

#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooRealProxy.h"

class RooSkewNormal : public RooAbsPdf {
  public:
    RooSkewNormal():RooAbsPdf(){}
    RooSkewNormal(const char *name, const char *title,
        RooAbsReal& _x, RooAbsReal& _location, RooAbsReal& _scale,
        RooAbsReal& _skew);
    RooSkewNormal(const RooSkewNormal& other, const char *name=0);
    virtual TObject* clone(const char* newname) const {
      return new RooSkewNormal(*this, newname);
    }
    virtual ~RooSkewNormal() {}

    /// Get value of PDF at current parameters (non--normalized)
    Double_t evaluate() const;

    /// Identify supported analytical integrals
    Int_t getAnalyticalIntegral(
        RooArgSet& allVars, RooArgSet& analVars, const char *rangeName) const;

    /// Compute analytical integrals
    Double_t analyticalIntegral(Int_t code, const char *rangeName) const;

  private:
    RooRealProxy x;
    RooRealProxy location;
    RooRealProxy scale;
    RooRealProxy skew;
    // generate CInt dictionaries
    ClassDef(RooSkewNormal,1);
};

#endif /* end of include guard: TauAnalysis_FittingTools_RooSkewNormal_h */
