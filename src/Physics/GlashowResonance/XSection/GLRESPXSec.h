//____________________________________________________________________________
/*!

\class    genie::GLRESPXSec

\brief    Nuebar cross section at the Glashow resonance (nuebar + e- -> W-).
          Is a concrete implementation of the XSecAlgorithmI interface.

\ref      T.K.Gaisser, F.Halzen and T.Stanev, Physics Reports 258:173 (1995)

\author   Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
          University of Liverpool & STFC Rutherford Appleton Laboratory

\created  May 04, 2005

\cpright  Copyright (c) 2003-2020, The GENIE Collaboration
          For the full text of the license visit http://copyright.genie-mc.org          
*/
//____________________________________________________________________________

#ifndef _GLASHOW_RESONANCE_PXSEC_H_
#define _GLASHOW_RESONANCE_PXSEC_H_

#include "Framework/EventGen/XSecAlgorithmI.h"

namespace genie {

class XSecIntegratorI;

class GLRESPXSec : public XSecAlgorithmI {

public:
  GLRESPXSec ();
  GLRESPXSec (string config);
  virtual ~GLRESPXSec ();

  // XSecAlgorithmI interface implementation
  double XSec            (const Interaction * i, KinePhaseSpace_t k) const;
  double Integral        (const Interaction * i) const;
  bool   ValidProcess    (const Interaction * i) const;

  // overload the Algorithm::Configure() methods to load private data
  // members from configuration options
  void Configure(const Registry & config);
  void Configure(string config);

private:
  void   LoadConfig (void);

  const XSecIntegratorI *        fXSecIntegrator;     ///< diff. xsec integrator

  double fWmin;            ///< Minimum value of W

};

}       // genie namespace

#endif  // _GLASHOW_RESONANCE_XSEC_H_
