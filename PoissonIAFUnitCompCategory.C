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
#include "PoissonIAFUnitCompCategory.h"
#include "NDPairList.h"
#include "CG_PoissonIAFUnitCompCategory.h"

PoissonIAFUnitCompCategory::PoissonIAFUnitCompCategory(Simulation& sim, const std::string& modelName, const NDPairList& ndpList) 
  : CG_PoissonIAFUnitCompCategory(sim, modelName, ndpList)
{
}
