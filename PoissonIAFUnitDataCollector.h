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

#ifndef PoissonIAFUnitDataCollector_H
#define PoissonIAFUnitDataCollector_H

#include "Lens.h"
#include "CG_PoissonIAFUnitDataCollector.h"
#include <memory>
#include <fstream>
#include <iostream>

class PoissonIAFUnitDataCollector : public CG_PoissonIAFUnitDataCollector
{
 public:
  void initialize(RNG& rng);
  void finalize(RNG& rng);
  virtual void dataCollection(Trigger* trigger, NDPairList* ndPairList);
  virtual void getNodeIndices(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_PoissonIAFUnitDataCollectorInAttrPSet* CG_inAttrPset, CG_PoissonIAFUnitDataCollectorOutAttrPSet* CG_outAttrPset);
  PoissonIAFUnitDataCollector();
  virtual ~PoissonIAFUnitDataCollector();
  virtual void duplicate(std::auto_ptr<PoissonIAFUnitDataCollector>& dup) const;
  virtual void duplicate(std::auto_ptr<Variable>& dup) const;
  virtual void duplicate(std::auto_ptr<CG_PoissonIAFUnitDataCollector>& dup) const;
 private:
  std::ofstream* spikes_file;
};

#endif
