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

#ifndef MihalasNieburSynapseIAFUnitCompCategory_H
#define MihalasNieburSynapseIAFUnitCompCategory_H

#include "Lens.h"
#include "CG_MihalasNieburSynapseIAFUnitCompCategory.h"

class NDPairList;

class MihalasNieburSynapseIAFUnitCompCategory : public CG_MihalasNieburSynapseIAFUnitCompCategory
{
 public:
  MihalasNieburSynapseIAFUnitCompCategory(Simulation& sim, const std::string& modelName, const NDPairList& ndpList);
  void initializeShared(RNG& rng);
 private:
  std::ofstream* AMPAindexs_file;
  std::ofstream* NMDARindexs_file;
  std::ofstream* GABARindexs_file;
  std::ostringstream os_AMPAindexs;
  std::ostringstream os_NMDARindexs;
  std::ostringstream os_GABARindexs;
};

#endif
