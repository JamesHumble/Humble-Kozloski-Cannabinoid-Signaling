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
#include "MihalasNieburSynapseIAFUnitCompCategory.h"
#include "NDPairList.h"
#include "CG_MihalasNieburSynapseIAFUnitCompCategory.h"
#include <fstream>
#include <sstream>

#define SHD getSharedMembers()

MihalasNieburSynapseIAFUnitCompCategory::MihalasNieburSynapseIAFUnitCompCategory(Simulation& sim, const std::string& modelName, const NDPairList& ndpList) 
  : CG_MihalasNieburSynapseIAFUnitCompCategory(sim, modelName, ndpList)
{
}

void MihalasNieburSynapseIAFUnitCompCategory::initializeShared(RNG& rng)
{
  std::ostringstream sysCall;
  sysCall<<"mkdir -p "<<SHD.sharedDirectory.c_str()<<";";
  try {
    int systemRet = system(sysCall.str().c_str());
    if (systemRet == -1)
      throw;
  } catch(...) {};
  if (SHD.op_saveIndexs)
    {
      int rank=getSimulation().getRank();
      int n=0;
      // Take it in turn opening and creating the file to create the stream on each node
      while (n<getSimulation().getNumProcesses()) {
        if (n==rank) {
          // AMPA
          os_AMPAindexs<<SHD.sharedDirectory<<"OutputAMPAIndexs"<<SHD.sharedFileExt;
          AMPAindexs_file=new std::ofstream(os_AMPAindexs.str().c_str(),
                                            std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
          AMPAindexs_file->close();
          // NMDAR
          os_NMDARindexs<<SHD.sharedDirectory<<"OutputNMDARIndexs"<<SHD.sharedFileExt;
          NMDARindexs_file=new std::ofstream(os_NMDARindexs.str().c_str(),
                                            std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
          NMDARindexs_file->close();
          // GABAR
          os_GABARindexs<<SHD.sharedDirectory<<"OutputGABARIndexs"<<SHD.sharedFileExt;
          GABARindexs_file=new std::ofstream(os_GABARindexs.str().c_str(),
                                             std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
          GABARindexs_file->close();
        }
        ++n;
        MPI_Barrier(MPI_COMM_WORLD); // wait node creating the stream to finish
      }
      // Now take it in turn writing to the file
      n=0;
      while (n<getSimulation().getNumProcesses()) {
        if (n==rank) {
          ShallowArray<MihalasNieburSynapseIAFUnit>::iterator it = _nodes.begin();
          ShallowArray<MihalasNieburSynapseIAFUnit>::iterator end = _nodes.end();
          AMPAindexs_file->open(os_AMPAindexs.str().c_str(),
                                std::ofstream::out | std::ofstream::app | std::ofstream::binary);
          NMDARindexs_file->open(os_NMDARindexs.str().c_str(),
                                std::ofstream::out | std::ofstream::app | std::ofstream::binary);
          GABARindexs_file->open(os_GABARindexs.str().c_str(),
                                 std::ofstream::out | std::ofstream::app | std::ofstream::binary);
          for (; it != end; ++it)
            {
              it->outputAMPAIndexs(*(AMPAindexs_file));
              it->outputNMDARIndexs(*(NMDARindexs_file));
              it->outputGABARIndexs(*(GABARindexs_file));
            }
          AMPAindexs_file->close();
          NMDARindexs_file->close();
          GABARindexs_file->close();
        }
        ++n;
        MPI_Barrier(MPI_COMM_WORLD); // wait for node writing to finish
      }
    }
}
