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
#include "CleftAstrocyteIAFUnit.h"
#include "CG_CleftAstrocyteIAFUnit.h"
#include "rndm.h"

#define SHD getSharedMembers()

void CleftAstrocyteIAFUnit::initialize(RNG& rng)
{
  // Check if not connected
  if (SHD.op_glutamateInput_check && !glutamateInput_check)
    std::cerr << "CleftAstrocyteIAFUnit: glutamateInput is not connected (possible errors)." << std::endl;
  if (SHD.op_GABAInput_check && !GABAInput_check)
    std::cerr << "CleftAstrocyteIAFUnit: GABAInput is not connected (possible errors)." << std::endl;
  if (SHD.op_glutamateeCBInput_check && !glutamateeCBInput_check)
    std::cerr << "CleftAstrocyteIAFUnit: glutamateeCBInput is not connected (possible errors)." << std::endl;
  if (SHD.op_GABAeCBInput_check && !GABAeCBInput_check)
    std::cerr << "CleftAstrocyteIAFUnit: GABAeCBInput is not connected (possible errors)." << std::endl;
  // Default starting values
  glutamate = 0.0;
  GABA = 0.0;
}

void CleftAstrocyteIAFUnit::update(RNG& rng)
{
  // Increase glutamate concentration in the cleft due to pre-synaptic release, and ...
  // ... decrease glutamate concentration due to astrocyte reuptake with GLT-1
  if (glutamateInput_check)
    glutamate += ((-glutamate + *(neurotransmitterInput.glutamate))
                  / SHD.glutamateDecayTau) * SHD.deltaT;

  // Increase GABA concentration in the cleft due to pre-synaptic release, and ...
  // ... decrease GABA concentration due to astrocyte reuptake with GLT-1
  if (GABAInput_check)
    GABA += ((-GABA + *(neurotransmitterInput.GABA))
             / SHD.GABADecayTau) * SHD.deltaT;

  // Combined spine and dendrite eCB into eCB concentrations bound for the glutamatergic and GABAergic boutons
  if (glutamateeCBInput_check && !GABAeCBInput_check)
    glutamateeCB = *(glutamateeCBInput.eCB);
  else if (glutamateeCBInput_check && GABAeCBInput_check)
    glutamateeCB = (*(glutamateeCBInput.eCB) * (1.0 - SHD.eCBCrossTalk))  
                     + (*(GABAeCBInput.eCB) * SHD.eCBCrossTalk);
  
  if (GABAeCBInput_check && !glutamateeCBInput_check)
    GABAeCB = *(GABAeCBInput.eCB);
  else if (GABAeCBInput_check && glutamateeCBInput_check)
    GABAeCB = (*(GABAeCBInput.eCB) * (1.0 - SHD.eCBCrossTalk))
      + (*(glutamateeCBInput.eCB) * SHD.eCBCrossTalk);  
}

void CleftAstrocyteIAFUnit::check_glutamateInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_CleftAstrocyteIAFUnitInAttrPSet* CG_inAttrPset, CG_CleftAstrocyteIAFUnitOutAttrPSet* CG_outAttrPset)
{
  if (SHD.op_glutamateInput_check && glutamateInput_check)
    std::cerr << "CleftAstrocyteIAFUnit: glutamateInput is already connected." << std::endl;
  glutamateInput_check = 1;  
}

void CleftAstrocyteIAFUnit::check_GABAInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_CleftAstrocyteIAFUnitInAttrPSet* CG_inAttrPset, CG_CleftAstrocyteIAFUnitOutAttrPSet* CG_outAttrPset)
{
  if (SHD.op_GABAInput_check && GABAInput_check)
    std::cerr << "CleftAstrocyteIAFUnit: GABAInput is already connected." << std::endl;
  GABAInput_check = 1;  
}

void CleftAstrocyteIAFUnit::check_glutamateeCBInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_CleftAstrocyteIAFUnitInAttrPSet* CG_inAttrPset, CG_CleftAstrocyteIAFUnitOutAttrPSet* CG_outAttrPset)
{
  if (SHD.op_glutamateeCBInput_check && glutamateeCBInput_check)
    std::cerr << "CleftAstrocyteIAFUnit: glutamateeCBInput is already connected." << std::endl;
  glutamateeCBInput_check = 1;  
}

void CleftAstrocyteIAFUnit::check_GABAeCBInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_CleftAstrocyteIAFUnitInAttrPSet* CG_inAttrPset, CG_CleftAstrocyteIAFUnitOutAttrPSet* CG_outAttrPset)
{
  if (SHD.op_GABAeCBInput_check && GABAeCBInput_check)
    std::cerr << "CleftAstrocyteIAFUnit: GABAeCBInput is already connected." << std::endl;
  GABAeCBInput_check = 1;  
}

CleftAstrocyteIAFUnit::~CleftAstrocyteIAFUnit()
{
}
