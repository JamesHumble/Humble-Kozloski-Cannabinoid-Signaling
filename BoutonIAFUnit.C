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
#include "BoutonIAFUnit.h"
#include "CG_BoutonIAFUnit.h"
#include "rndm.h"
#include <fstream>
#include <sstream>

#define SHD getSharedMembers()

void BoutonIAFUnit::initialize(RNG& rng)
{
  // Check if not connected
  if (SHD.op_spikeInput_check && !spikeInput_check)
    std::cerr << "BoutonIAFUnit: spikeInput is not connected (possible errors)." << std::endl;
  if (SHD.op_eCBInput_check && !eCBInput_check)
    std::cerr << "BoutonIAFUnit: eCBInput is not connected (possible errors)." << std::endl;
  if (SHD.op_CB1RInput_check && !CB1RInput_check)
    std::cerr << "BoutonIAFUnit: CB1RInput is not connected (possible errors)." << std::endl;
  // Default starting values
  neurotransmitter = 0.0;
  CB1R = 0.0;
  CB1Runbound = 0.0;  
  CB1Rrise = 0.0;
  CB1Rcurrent = 0.0;
  MAGLrise = 0.0;
  MAGL = 0.0;
}

void BoutonIAFUnit::update(RNG& rng)
{
  // ##### Neurotransmitter release #####
  // If there is a spike, release neurotransmitter
  neurotransmitter = *(spikeInput.spike) ? availableNeurotransmitter : 0.0;



  // ##### CB1R #####
  CB1R = *(CB1RInput.Y) * SHD.CB1RInputWeight; // weight is scaling from Goodwin model Y      
  if (CB1R > 1.0) // bound just in case
    CB1R = 1.0;
  else if (CB1R < 0.0)
    CB1R = 0.0;

  CB1Runbound = CB1R - *(eCBInput.eCB);
  if (CB1Runbound < 0.0)
    CB1Runbound = 0.0;
  if (neurotransmitterType == 1 || SHD.op_saveCB1Rcurrent) // Only for GABA or if glutamate and saving
    {
      CB1Rrise += ((-CB1Rrise + std::min(*(eCBInput.eCB), CB1R)) / SHD.CB1RriseTau) * SHD.deltaT;
      CB1Rcurrent += ((-CB1Rcurrent + CB1Rrise) / SHD.CB1RfallTau) * SHD.deltaT;
    }
  


  // ##### Inhibit neurotransmitter #####
  // Recovery neurotransmitter
  availableNeurotransmitter += ((maxNeurotransmitter - availableNeurotransmitter)
                                / SHD.neurotransmitterRecoverTau[neurotransmitterType]) * SHD.deltaT;
  // Inhibit the neurotransmitter release with the quantity of eCB and CB1R, i.e. the current amount bound only, minimum
  availableNeurotransmitter -= SHD.neurotransmitterAdaptRate[neurotransmitterType]
    * std::min(*(eCBInput.eCB), CB1R);
  // Limit neurotransmitter to >= 0
  if (availableNeurotransmitter < 0.0)
    availableNeurotransmitter = 0.0;



  // ##### MAGL #####
  if (neurotransmitterType == 1) // Only for GABAergic boutons
    {
      MAGLrise += ((-MAGLrise + (CB1Rcurrent * SHD.MAGLsensitivity))
                   / SHD.MAGLriseTau) * SHD.deltaT;
      MAGL += ((-MAGL + MAGLrise) / SHD.MAGLfallTau) * SHD.deltaT;      
    }
}

void BoutonIAFUnit::outputIndexs(std::ofstream& fs)
{
  fs.write(reinterpret_cast<char *>(&(spikeInput.col)), sizeof(spikeInput.col));
  fs.write(reinterpret_cast<char *>(&(spikeInput.row)), sizeof(spikeInput.row));
}

void BoutonIAFUnit::setSpikeIndices(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_BoutonIAFUnitInAttrPSet* CG_inAttrPset, CG_BoutonIAFUnitOutAttrPSet* CG_outAttrPset)
{
  spikeInput.row =  getIndex()+1; // +1 is for Matlab
  spikeInput.col = CG_node->getIndex()+1;
}

void BoutonIAFUnit::check_spikeInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_BoutonIAFUnitInAttrPSet* CG_inAttrPset, CG_BoutonIAFUnitOutAttrPSet* CG_outAttrPset)
{
  if (SHD.op_spikeInput_check && spikeInput_check)
    std::cerr << "BoutonIAFUnit: spikeInput is already connected." << std::endl;
  spikeInput_check = 1;
}

void BoutonIAFUnit::check_eCBInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_BoutonIAFUnitInAttrPSet* CG_inAttrPset, CG_BoutonIAFUnitOutAttrPSet* CG_outAttrPset)
{
  if (SHD.op_eCBInput_check && eCBInput_check)
    std::cerr << "BoutonIAFUnit: eCBInput is already connected." << std::endl;
  eCBInput_check = 1;
}

void BoutonIAFUnit::check_CB1RInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_BoutonIAFUnitInAttrPSet* CG_inAttrPset, CG_BoutonIAFUnitOutAttrPSet* CG_outAttrPset)
{
  if (SHD.op_CB1RInput_check && CB1RInput_check)
    std::cerr << "BoutonIAFUnit: CB1RInput is already connected." << std::endl;
  CB1RInput_check = 1;
}

BoutonIAFUnit::~BoutonIAFUnit()
{
}
