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
#include "DendriteIAFUnit.h"
#include "CG_DendriteIAFUnit.h"
#include "rndm.h"
#include <fstream>
#include <sstream>

#define SHD getSharedMembers()
#define ITER getSimulation().getIteration()

void DendriteIAFUnit::initialize(RNG& rng)
{
  // Check if not connected
  if (SHD.op_neurotransmitterInput_check && !neurotransmitterInput_check)
    std::cerr << "DendriteIAFUnit: neurotransmitterInputg is not connected (possible errors)." << std::endl;
  if (SHD.op_postSpikeInput_check && !postSpikeInput_check)
    std::cerr << "DendriteIAFUnit: postSpikeInputg is not connected (possible errors)." << std::endl;
  // Default starting values
  GABArise = 0.0;
  GABAcurrent = 0.0;
  Carise = 0.0;
  Ca = 0.0;
  eCB = 0.0;
}

void DendriteIAFUnit::update(RNG& rng)
{  
  // If the simulation has reached a certain period, apply a perturbation
  if (SHD.op_perturbation && ITER == SHD.perturbationT)
    GABAweight = drandom(0.0, 1.5, rng);
  


  // ##### GABA #####
  // GABA input is the minimum of GABA weight and neurotransmitter
  // Only updated when there is a pre-spike
  double GABAinput = std::min(*(neurotransmitterInput.neurotransmitter), GABAweight);

  // Update GABA rise with the GABA activity
  GABArise += ((-GABArise + GABAinput) / SHD.GABAriseTau ) * SHD.deltaT;
  // Update GABA current (fall) with the GABA rise
  GABAcurrent += ((-GABAcurrent + GABArise) / SHD.GABAfallTau) * SHD.deltaT;



  // ##### Ca2+ #####
  // Update Ca2+ rise with VSCC and BP
  Carise += ((-Carise + (*(postSpikeInput.spike) * SHD.CaVSCC))
             / SHD.CariseTau) * SHD.deltaT;
  // Update Ca2+ fall with Ca2+ rise
  Ca += ((-Ca + Carise) / SHD.CafallTau) * SHD.deltaT;



  // ##### endocannabinoids #####
  // Update eCB (is always in the range 0 to 1)
  // with Ca2+
  eCB = eCBproduction(Ca);
}

void DendriteIAFUnit::outputWeights(std::ofstream& fs)
{
  float temp = (float) GABAweight;
  fs.write(reinterpret_cast<char *>(&temp), sizeof(temp));
}

void DendriteIAFUnit::check_neurotransmitterInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_DendriteIAFUnitInAttrPSet* CG_inAttrPset, CG_DendriteIAFUnitOutAttrPSet* CG_outAttrPset)
{
  if (SHD.op_neurotransmitterInput_check && neurotransmitterInput_check)
    std::cerr << "DendriteIAFUnit: neurotransmitterInput is already connected." << std::endl;
  neurotransmitterInput_check = 1;
}

void DendriteIAFUnit::check_postSpikeInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_DendriteIAFUnitInAttrPSet* CG_inAttrPset, CG_DendriteIAFUnitOutAttrPSet* CG_outAttrPset)
{
  if (SHD.op_postSpikeInput_check && postSpikeInput_check)
    std::cerr << "DendriteIAFUnit: postSpikeInput is already connected." << std::endl;
  postSpikeInput_check = 1;
}

DendriteIAFUnit::~DendriteIAFUnit()
{
}

double DendriteIAFUnit::eCBsigmoid(double Ca)
{
  return 1.0 / ( 1.0 + exp(-SHD.eCBprodC * (Ca - SHD.eCBprodD)) );
}

double DendriteIAFUnit::eCBproduction(double Ca)
{
  // Computes the sigmoidal production of cannabinoids depending on Ca2+
  // NOTE: this is mirrored in DendriteIAFUnitDataCollector. If changed here, change there too.
  double eCB = 0.0;
  // 1. the general sigmoid
  eCB = eCBsigmoid(Ca);
  // 2. make zero eCB at zero Ca2+
  eCB -= eCBsigmoid(0.0);
  // 3. Make one eCB at >= one Ca2+
  eCB *= 1.0 / (eCBsigmoid(1.0) - eCBsigmoid(0.0));
  if (eCB > 1.0)
    eCB = 1.0;

  return eCB;
}
