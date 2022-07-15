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
#include "SpineIAFUnit.h"
#include "CG_SpineIAFUnit.h"
#include "rndm.h"
#include <fstream>
#include <sstream>

#define SHD getSharedMembers()
#define ITER getSimulation().getIteration()

void SpineIAFUnit::initialize(RNG& rng)
{
  // Check if not connected
  if (SHD.op_neurotransmitterInput_bouton_check && !neurotransmitterInput_bouton_check)
    std::cerr << "SpineIAFUnit: neurotransmitterInput_bouton is not connected (possible errors)." << std::endl;
  if (SHD.op_neurotransmitterInput_cleft_check && !neurotransmitterInput_cleft_check)
    std::cerr << "SpineIAFUnit: neurotransmitterInput_cleft is not connected (possible errors)." << std::endl;
  if (SHD.op_postSpikeInput_check && !postSpikeInput_check)
    std::cerr << "SpineIAFUnit: postSpikeInput is not connected (possible errors)." << std::endl;
  // Default starting values
  AMPArise = 0.0;
  AMPAcurrent = 0.0;
  mGluR5rise = 0.0;
  mGluR5current = 0.0;
  NMDARopen = 0.0;
  NMDARrise = 0.0;
  NMDARcurrent = 0.0;
  Carise = 0.0;
  Ca = 0.0;
  eCB = 0.0;
}

void SpineIAFUnit::update(RNG& rng)
{  
  // If the simulation has reached a certain period, apply a perturbation
  if (SHD.op_perturbation && ITER == SHD.perturbationT)
    AMPAweight = drandom(0.0, 1.5, rng);
  


  // ##### AMPA #####
  // AMPA input is the minimum of AMPA weight and neurotransmitter from the bouton
  // Only updated when there is a pre-spike
  double AMPAinput = std::min(*(neurotransmitterInput_bouton.neurotransmitter), AMPAweight);

  // Update AMPA rise with the AMPA activity
  AMPArise += ((-AMPArise + AMPAinput) / SHD.AMPAriseTau ) * SHD.deltaT;
  // Update AMPA current (fall) with the AMPA rise
  AMPAcurrent += ((-AMPAcurrent + AMPArise) / SHD.AMPAfallTau) * SHD.deltaT;



  // ##### mGluR5 #####
  // mGluR5 input is any excess neurotransmitter from the bouton bigger than AMPA weight
  double mGluR5input = 0.0;
  if (*(neurotransmitterInput_bouton.neurotransmitter) > AMPAweight)
    mGluR5input = (*(neurotransmitterInput_bouton.neurotransmitter) - AMPAweight) * SHD.mGluR5sensitivity; // adjust the sensitivity as well

  // Update mGluR5 rise with the mGluR5 activity
  mGluR5rise += ((-mGluR5rise + mGluR5input) / SHD.mGluR5riseTau ) * SHD.deltaT;
  // Update mGluR5 current (fall) with the mGluR5 rise
  mGluR5current += ((-mGluR5current + mGluR5rise) / SHD.mGluR5fallTau) * SHD.deltaT;



  // ##### NMDAR #####
  double NMDARinput = *(neurotransmitterInput_cleft.neurotransmitter)
    - *(neurotransmitterInput_bouton.neurotransmitter);
  NMDARopen += ((-NMDARopen + (NMDARinput > 0.0 ? NMDARinput : 0.0)) // only what is unbound: N.B. glutamate binds and unbinds from AMPA and mGluR5 in one dt!
                / SHD.NMDARopenTau) * SHD.deltaT;
      
  // Update NMDAR rise with the NMDAR activity
  NMDARrise += ((-NMDARrise + (*(postSpikeInput.spike) * NMDARopen * NMDARweight))
                / SHD.NMDARriseTau ) * SHD.deltaT;
  // Update NMDAR current (fall) with the NMDAR rise
  NMDARcurrent += ((-NMDARcurrent + NMDARrise) / SHD.NMDARfallTau) * SHD.deltaT;


  
  // ##### Ca2+ #####
  // Ca2+ input
  double CaVSCCinput = 0.0;
  if (SHD.op_CaVSCCdepend)
    CaVSCCinput = SHD.CaVSCC * pow(AMPAweight, SHD.CaVSCCpow);
  else
    CaVSCCinput = SHD.CaVSCC;
  double CaInput = (CaVSCCinput * AMPAcurrent) + (SHD.NMDARCasensitivity * NMDARcurrent);

  // Update Ca2+ rise with VSCC and BP
  Carise += ((-Carise + CaInput) / SHD.CariseTau) * SHD.deltaT;
  // Update Ca2+ fall with Ca2+ rise
  Ca += ((-Ca + Carise) / SHD.CafallTau) * SHD.deltaT;



  // ##### endocannabinoids #####
  // Update eCB (is always in the range 0 to 1)
  // with Ca2+ and the mGluR5 modulation (AND gate)
  eCB = eCBproduction(Ca * mGluR5modulation(mGluR5current));
}

void SpineIAFUnit::outputAMPAweights(std::ofstream& fs)
{
  float temp = (float) AMPAweight;
  fs.write(reinterpret_cast<char *>(&temp), sizeof(temp));
}

void SpineIAFUnit::outputNMDARweights(std::ofstream& fs)
{
  float temp = (float) NMDARweight;
  fs.write(reinterpret_cast<char *>(&temp), sizeof(temp));
}

void SpineIAFUnit::check_neurotransmitterInput_bouton(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_SpineIAFUnitInAttrPSet* CG_inAttrPset, CG_SpineIAFUnitOutAttrPSet* CG_outAttrPset)
{  
  if (SHD.op_neurotransmitterInput_bouton_check && neurotransmitterInput_bouton_check)
    std::cerr << "SpineIAFUnit: neurotransmitterInput_bouton is already connected." << std::endl;
  neurotransmitterInput_bouton_check = 1;
}

void SpineIAFUnit::check_neurotransmitterInput_cleft(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_SpineIAFUnitInAttrPSet* CG_inAttrPset, CG_SpineIAFUnitOutAttrPSet* CG_outAttrPset)
{  
  if (SHD.op_neurotransmitterInput_cleft_check && neurotransmitterInput_cleft_check)
    std::cerr << "SpineIAFUnit: neurotransmitterInput_cleft is already connected." << std::endl;
  neurotransmitterInput_cleft_check = 1;  
}

void SpineIAFUnit::check_postSpikeInput(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_SpineIAFUnitInAttrPSet* CG_inAttrPset, CG_SpineIAFUnitOutAttrPSet* CG_outAttrPset)
{
  if (SHD.op_postSpikeInput_check && postSpikeInput_check)
    std::cerr << "SpineIAFUnit: postSpikeInput is already connected." << std::endl;
  postSpikeInput_check = 1;
}

SpineIAFUnit::~SpineIAFUnit()
{
}

double SpineIAFUnit::eCBsigmoid(double Ca)
{
  return 1.0 / ( 1.0 + exp(-SHD.eCBprodC * (Ca - SHD.eCBprodD)) );
}

double SpineIAFUnit::eCBproduction(double Ca)
{
  // Computes the sigmoidal production of cannabinoids depending on Ca2+
  // NOTE: this is mirrored in SpineIAFUnitDataCollector. If changed here, change there too.
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

double SpineIAFUnit::mGluR5modulation(double mGluR5)
{
  return eCBproduction(mGluR5); // just use the same modified sigmoid
}
