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

#ifndef DendriteIAFUnitDataCollector_H
#define DendriteIAFUnitDataCollector_H

#include "Lens.h"
#include "CG_DendriteIAFUnitDataCollector.h"
#include <memory>
#include <fstream>
#include <iostream>

class DendriteIAFUnitDataCollector : public CG_DendriteIAFUnitDataCollector
{
 public:
  void initialize(RNG& rng);
  void finalize(RNG& rng);
  virtual void dataCollection(Trigger* trigger, NDPairList* ndPairList);
  virtual void getNodeIndices(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_DendriteIAFUnitDataCollectorInAttrPSet* CG_inAttrPset, CG_DendriteIAFUnitDataCollectorOutAttrPSet* CG_outAttrPset);
  DendriteIAFUnitDataCollector();
  virtual ~DendriteIAFUnitDataCollector();
  virtual void duplicate(std::auto_ptr<DendriteIAFUnitDataCollector>& dup) const;
  virtual void duplicate(std::auto_ptr<Variable>& dup) const;
  virtual void duplicate(std::auto_ptr<CG_DendriteIAFUnitDataCollector>& dup) const;
 private:
  std::ofstream* GABA_file;
  std::ofstream* Ca_file;
  std::ofstream* eCBproduction_file; // just for the production function
  std::ofstream* eCB_file;
  double eCBsigmoid(double Ca);
  double eCBproduction(double Ca);
};

#endif
