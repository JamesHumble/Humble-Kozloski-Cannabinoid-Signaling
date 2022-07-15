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

#ifndef SpineIAFUnit_H
#define SpineIAFUnit_H

#include "Lens.h"
#include "CG_SpineIAFUnit.h"
#include "rndm.h"

class SpineIAFUnit : public CG_SpineIAFUnit
{
 public:
  void initialize(RNG& rng);
  void update(RNG& rng);
  void outputAMPAweights(std::ofstream& fs);
  void outputNMDARweights(std::ofstream& fs);
  virtual void check_neurotransmitterInput_bouton(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_SpineIAFUnitInAttrPSet* CG_inAttrPset, CG_SpineIAFUnitOutAttrPSet* CG_outAttrPset);
  virtual void check_neurotransmitterInput_cleft(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_SpineIAFUnitInAttrPSet* CG_inAttrPset, CG_SpineIAFUnitOutAttrPSet* CG_outAttrPset);
  virtual void check_postSpikeInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_SpineIAFUnitInAttrPSet* CG_inAttrPset, CG_SpineIAFUnitOutAttrPSet* CG_outAttrPset);
  virtual ~SpineIAFUnit();
 private:
  double eCBsigmoid(double Ca);
  double eCBproduction(double Ca);
  double mGluR5modulation(double mGluR5);
};

#endif
