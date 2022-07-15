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

#ifndef CleftAstrocyteIAFUnitCompCategory_H
#define CleftAstrocyteIAFUnitCompCategory_H

#include "Lens.h"
#include "CG_CleftAstrocyteIAFUnitCompCategory.h"

class NDPairList;

class CleftAstrocyteIAFUnitCompCategory : public CG_CleftAstrocyteIAFUnitCompCategory
{
 public:
  CleftAstrocyteIAFUnitCompCategory(Simulation& sim, const std::string& modelName, const NDPairList& ndpList);
  void initializeShared(RNG& rng);
};

#endif
