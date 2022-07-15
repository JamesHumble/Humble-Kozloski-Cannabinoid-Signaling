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

#ifndef DendriteIAFUnitCompCategory_H
#define DendriteIAFUnitCompCategory_H

#include "Lens.h"
#include "CG_DendriteIAFUnitCompCategory.h"

class NDPairList;

class DendriteIAFUnitCompCategory : public CG_DendriteIAFUnitCompCategory
{
 public:
  DendriteIAFUnitCompCategory(Simulation& sim, const std::string& modelName, const NDPairList& ndpList);
  void initializeShared(RNG& rng);
  void outputWeightsShared(RNG& rng);
 private:
  std::ofstream* weights_file;
  std::ostringstream os_weights;
};

#endif
