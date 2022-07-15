// =================================================================
// Licensed Materials - Property of IBM
//
// "Restricted Materials of IBM"
//
// BCM-YKT-07-18-2017
//
// (C) Copyright IBM Corp. 2005-2022  All rights reserved
//
// US Government Users Restricted Rights -
// Use, duplication or disclosure restricted by
// GSA ADP Schedule Contract with IBM Corp.
//
// =================================================================

#ifndef BoutonIAFUnitCompCategory_H
#define BoutonIAFUnitCompCategory_H

#include "Lens.h"
#include "CG_BoutonIAFUnitCompCategory.h"

class NDPairList;

class BoutonIAFUnitCompCategory : public CG_BoutonIAFUnitCompCategory
{
 public:
  BoutonIAFUnitCompCategory(Simulation& sim, const std::string& modelName, const NDPairList& ndpList);
  void initializeShared(RNG& rng);
 private:
  std::ofstream** indexs_file;
  std::ostringstream* os_indexs;
};
#endif
