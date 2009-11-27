#ifndef TauAnalysis_FittingTools_TF1landauXgausWrapper_h
#define TauAnalysis_FittingTools_TF1landauXgausWrapper_h

/** \class TF1landauXgausWrapper
 *
 * Numerical implementation of Landau density convoluted with Gaussian
 * 
 * \author Christian Veelken, UC Davis
 *
 * \version $Revision: 1.2 $
 *
 * $Id: TF1landauXgausWrapper.h,v 1.2 2009/09/04 16:18:23 veelken Exp $
 *
 */

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/FittingTools/interface/TF1WrapperBase.h"

class TF1landauXgausWrapper : public TF1WrapperBase
{
 public:
  // constructor 
  explicit TF1landauXgausWrapper(const edm::ParameterSet& cfg);   
  
  // destructor
  ~TF1landauXgausWrapper() {}
};

#endif  

