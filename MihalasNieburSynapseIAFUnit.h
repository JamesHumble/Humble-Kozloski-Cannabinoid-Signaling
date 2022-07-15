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

#ifndef MihalasNieburSynapseIAFUnit_H
#define MihalasNieburSynapseIAFUnit_H

#include "Lens.h"
#include "CG_MihalasNieburSynapseIAFUnit.h"
#include "rndm.h"

class MihalasNieburSynapseIAFUnit : public CG_MihalasNieburSynapseIAFUnit
{
 public:
  void initialize(RNG& rng);
  void update(RNG& rng);
  void threshold(RNG& rng);
  void outputAMPAIndexs(std::ofstream& fs);
  void outputNMDARIndexs(std::ofstream& fs);
  void outputGABARIndexs(std::ofstream& fs);
  virtual void setAMPAIndices(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_MihalasNieburSynapseIAFUnitInAttrPSet* CG_inAttrPset, CG_MihalasNieburSynapseIAFUnitOutAttrPSet* CG_outAttrPset);
  virtual void setNMDARIndices(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_MihalasNieburSynapseIAFUnitInAttrPSet* CG_inAttrPset, CG_MihalasNieburSynapseIAFUnitOutAttrPSet* CG_outAttrPset);
  virtual void setGABARIndices(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_MihalasNieburSynapseIAFUnitInAttrPSet* CG_inAttrPset, CG_MihalasNieburSynapseIAFUnitOutAttrPSet* CG_outAttrPset);
  virtual ~MihalasNieburSynapseIAFUnit();
};

#endif
