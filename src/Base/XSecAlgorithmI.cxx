//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - May 03, 2004

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Base/XSecAlgorithmI.h"
#include "Messenger/Messenger.h"

using namespace genie;

//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI() :
Algorithm()
{

}
//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI(string name) :
Algorithm(name)
{

}
//___________________________________________________________________________
XSecAlgorithmI::XSecAlgorithmI(string name, string config) :
Algorithm(name, config)
{

}
//___________________________________________________________________________
XSecAlgorithmI::~XSecAlgorithmI()
{

}
//___________________________________________________________________________
bool XSecAlgorithmI::ValidKinematics(const Interaction* interaction) const
{
// can offer common implementation for all concrete x-section models because
// the input interaction is aware of its kinematic limits

  if(interaction->TestBit(kISkipKinematicChk)) return true;

  KPhaseSpace kps = interaction->PhaseSpace();

  if(!kps.IsAboveThreshold()) {
     LOG("XSecBase", pINFO)  << "*** Below energy threshold";
     return false;
  }
  if(!kps.IsAllowed()) {
     LOG("XSecBase", pINFO)  << "*** Not in allowed kinematical space";
     return false;
  }
  return true;
}
//___________________________________________________________________________

