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

#ifndef DendriteIAFUnit_H
#define DendriteIAFUnit_H

#include "Lens.h"
#include "CG_DendriteIAFUnit.h"
#include "rndm.h"

class DendriteIAFUnit : public CG_DendriteIAFUnit
{
 public:
  void initialize(RNG& rng);
  void update(RNG& rng);
  void outputWeights(std::ofstream& fs);
  virtual void check_neurotransmitterInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_DendriteIAFUnitInAttrPSet* CG_inAttrPset, CG_DendriteIAFUnitOutAttrPSet* CG_outAttrPset);
  virtual void check_postSpikeInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_DendriteIAFUnitInAttrPSet* CG_inAttrPset, CG_DendriteIAFUnitOutAttrPSet* CG_outAttrPset);
  virtual ~DendriteIAFUnit();
 private:
  double eCBsigmoid(double Ca);
  double eCBproduction(double Ca);
};

#endif
