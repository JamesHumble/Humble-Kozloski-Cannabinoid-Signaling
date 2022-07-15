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
#include "SpineIAFUnitCompCategory.h"
#include "NDPairList.h"
#include "CG_SpineIAFUnitCompCategory.h"
#include <fstream>
#include <sstream>

#define SHD getSharedMembers()
#define ITER getSimulation().getIteration()

SpineIAFUnitCompCategory::SpineIAFUnitCompCategory(Simulation& sim, const std::string& modelName, const NDPairList& ndpList)
  : CG_SpineIAFUnitCompCategory(sim, modelName, ndpList)
{
  SHD.collectAMPAweightsNext = 0;
  SHD.collectNMDARweightsNext = 0;
}

void SpineIAFUnitCompCategory::initializeShared(RNG& rng)
{
  // AMPA
  std::ostringstream sysCall;
  sysCall<<"mkdir -p "<<SHD.sharedDirectory.c_str()<<";";
  try {
    int systemRet = system(sysCall.str().c_str());
    if (systemRet == -1)
      throw;
  } catch(...) {};  
  if (SHD.op_saveAMPAweights)
    {
      int n=SHD.collectAMPAweightsOn.size();
      if (n>0)
        {
          int rank=getSimulation().getRank();          
          for (int i=0; i<n; i++)
            {
              int r=0;
              while (r<getSimulation().getNumProcesses())
                {                  
                  if (r==rank)
                    {
                      os_AMPAweights.str(std::string());
                      os_AMPAweights<<SHD.sharedDirectory<<SHD.sharedFilePrep
                                <<"SpineAMPAweights_"<<SHD.collectAMPAweightsOn[i]
                                <<SHD.sharedFileApp<<SHD.sharedFileExt;
                      AMPAweights_file=new std::ofstream(os_AMPAweights.str().c_str(),
                                                     std::ofstream::out | std::ofstream::trunc
                                                     | std::ofstream::binary);
                      AMPAweights_file->close();
                    }                  
                  ++r;
                  MPI_Barrier(MPI_COMM_WORLD); // wait node creating the stream to finish
                }
            }
        }
    }
  // NMDAR
  sysCall<<"mkdir -p "<<SHD.sharedDirectory.c_str()<<";";
  try {
    int systemRet = system(sysCall.str().c_str());
    if (systemRet == -1)
      throw;
  } catch(...) {};  
  if (SHD.op_saveNMDARweights)
    {
      int n=SHD.collectNMDARweightsOn.size();
      if (n>0)
        {
          int rank=getSimulation().getRank();          
          for (int i=0; i<n; i++)
            {
              int r=0;
              while (r<getSimulation().getNumProcesses())
                {                  
                  if (r==rank)
                    {
                      os_NMDARweights.str(std::string());
                      os_NMDARweights<<SHD.sharedDirectory<<SHD.sharedFilePrep
                                <<"SpineNMDARweights_"<<SHD.collectNMDARweightsOn[i]
                                <<SHD.sharedFileApp<<SHD.sharedFileExt;
                      NMDARweights_file=new std::ofstream(os_NMDARweights.str().c_str(),
                                                     std::ofstream::out | std::ofstream::trunc
                                                     | std::ofstream::binary);
                      NMDARweights_file->close();
                    }                  
                  ++r;
                  MPI_Barrier(MPI_COMM_WORLD); // wait node creating the stream to finish
                }
            }
        }
    }
}

void SpineIAFUnitCompCategory::outputAMPAweightsShared(RNG& rng)
{
  if (SHD.op_saveAMPAweights)
    {
      int n=SHD.collectAMPAweightsOn.size();
      if (SHD.collectAMPAweightsOn[SHD.collectAMPAweightsNext]==ITER)
        {
          os_AMPAweights.str(std::string());
          os_AMPAweights<<SHD.sharedDirectory<<SHD.sharedFilePrep
                    <<"SpineAMPAweights_"<<SHD.collectAMPAweightsOn[SHD.collectAMPAweightsNext]
                    <<SHD.sharedFileApp<<SHD.sharedFileExt;
          if (SHD.collectAMPAweightsOn.size()-1 > SHD.collectAMPAweightsNext)
            SHD.collectAMPAweightsNext++;
          int rank=getSimulation().getRank();
          int r=0;
          while (r<getSimulation().getNumProcesses())
            {
              if (r==rank) {
                ShallowArray<SpineIAFUnit>::iterator it = _nodes.begin();
                ShallowArray<SpineIAFUnit>::iterator end = _nodes.end();
                AMPAweights_file->open(os_AMPAweights.str().c_str(),
                                   std::ofstream::out | std::ofstream::app | std::ofstream::binary);
                for (; it != end; ++it)
                  (*it).outputAMPAweights(*AMPAweights_file);
                AMPAweights_file->close();
              }
              ++r;
              MPI_Barrier(MPI_COMM_WORLD); // wait for node writing to finish
            }
        }
    }
}

void SpineIAFUnitCompCategory::outputNMDARweightsShared(RNG& rng)
{
  if (SHD.op_saveNMDARweights)
    {
      int n=SHD.collectNMDARweightsOn.size();
      if (SHD.collectNMDARweightsOn[SHD.collectNMDARweightsNext]==ITER)
        {
          os_NMDARweights.str(std::string());
          os_NMDARweights<<SHD.sharedDirectory<<SHD.sharedFilePrep
                    <<"SpineNMDARweights_"<<SHD.collectNMDARweightsOn[SHD.collectNMDARweightsNext]
                    <<SHD.sharedFileApp<<SHD.sharedFileExt;
          if (SHD.collectNMDARweightsOn.size()-1 > SHD.collectNMDARweightsNext)
            SHD.collectNMDARweightsNext++;
          int rank=getSimulation().getRank();
          int r=0;
          while (r<getSimulation().getNumProcesses())
            {
              if (r==rank) {
                ShallowArray<SpineIAFUnit>::iterator it = _nodes.begin();
                ShallowArray<SpineIAFUnit>::iterator end = _nodes.end();
                NMDARweights_file->open(os_NMDARweights.str().c_str(),
                                   std::ofstream::out | std::ofstream::app | std::ofstream::binary);
                for (; it != end; ++it)
                  (*it).outputNMDARweights(*NMDARweights_file);
                NMDARweights_file->close();
              }
              ++r;
              MPI_Barrier(MPI_COMM_WORLD); // wait for node writing to finish
            }
        }
    }
}
