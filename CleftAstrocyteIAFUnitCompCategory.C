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

#include "Lens.h"
#include "CleftAstrocyteIAFUnitCompCategory.h"
#include "NDPairList.h"
#include "CG_CleftAstrocyteIAFUnitCompCategory.h"

#define SHD getSharedMembers()

CleftAstrocyteIAFUnitCompCategory::CleftAstrocyteIAFUnitCompCategory(Simulation& sim, const std::string& modelName, const NDPairList& ndpList) 
  : CG_CleftAstrocyteIAFUnitCompCategory(sim, modelName, ndpList)
{
}

void CleftAstrocyteIAFUnitCompCategory::initializeShared(RNG& rng) 
{
  if (SHD.eCBCrossTalk > 1.0 || SHD.eCBCrossTalk < 0.0)
    std::cerr << "CleftAstrocyteIAFUnitCompCategory: eCBCrossTalk is not a valid value." << std::endl;  
}
