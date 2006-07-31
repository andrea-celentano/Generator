//____________________________________________________________________________
/*!

\class   genie::PDGCodeList

\brief   A list of PDG codes

\author  Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory

\created May 13, 2005

\cpright  Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
          All rights reserved.
          For the licensing terms see $GENIE/USER_LICENSE.
*/
//____________________________________________________________________________

#ifndef _PDG_CODE_LIST_H_
#define _PDG_CODE_LIST_H_

#include <vector>
#include <ostream>

using std::vector;
using std::ostream;

namespace genie {

class PDGCodeList : public vector<int> {

public :

  PDGCodeList(bool allowdup=false);
  PDGCodeList(size_type n, bool allowdup=false);
  PDGCodeList(const PDGCodeList & list);
  ~PDGCodeList();

  //! override the vector<int> insertion methods to explicitly check for
  //! PDG code validity and that no PDG code is listed more than once
  void push_back  (int pdg_code);
  void insert     (iterator pos, size_type n, const int& x);

  //! PDG code checks used by PDGCodeList
  bool CheckPDGCode        (int pdg_code);
  bool ExistsInPDGLibrary  (int pdg_code);
  bool ExistsInPDGCodeList (int pdg_code);

  //! copy / print
  void Copy  (const PDGCodeList & list);
  void Print (ostream & stream) const;

  //! check state
  bool DuplEntriesAllowed(void) const { return fAllowDuplicateEntries; }

  //! overloaded operators
  PDGCodeList &    operator =  (const PDGCodeList & list);
  friend ostream & operator << (ostream & stream, const PDGCodeList & list);

private:

  bool fAllowDuplicateEntries; ///< allow duplicate entries in the list?
};

}      // genie namespace

#endif // _PDG_CODE_LIST_H_
