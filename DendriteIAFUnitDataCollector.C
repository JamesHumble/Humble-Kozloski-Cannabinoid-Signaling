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
#include "DendriteIAFUnitDataCollector.h"
#include "CG_DendriteIAFUnitDataCollector.h"
#include "NodeDescriptor.h"
#include "Node.h"
#include <map>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <utility>

void DendriteIAFUnitDataCollector::initialize(RNG& rng)
{
  // Sort pointers by indices, row major
  std::map<unsigned,
	   std::map<unsigned,
                    std::map<unsigned,
                             std::pair<
                               std::pair<double*, double*>, // second.first.first second.first.second
                               double* // second.second
                               >
                             >
                    >
           >
    sorter;
  assert(rows.size()==slices.size());
  assert(cols.size()==slices.size());
  assert(slices.size()==GABAcurrent.size());
  assert(slices.size()==Ca.size());
  assert(slices.size()==eCB.size());
  int sz=GABAcurrent.size();
  int mxrow=0;
  int mxcol=0;
  for (int j=0; j<sz; ++j)
    {
      sorter[rows[j]][cols[j]][slices[j]]=std::make_pair(
                                                         std::make_pair(GABAcurrent[j], Ca[j]),
                                                         eCB[j]
                                                         );
      if (mxrow<rows[j]) mxrow=rows[j];
      if (mxcol<cols[j]) mxcol=cols[j];
      if (mxslice<slices[j]) mxslice=slices[j];
    }
  GABAcurrent.clear();
  Ca.clear();
  eCB.clear();
  std::map<unsigned,
	   std::map<unsigned,
                    std::map<unsigned,
                             std::pair<
                               std::pair<double*, double*>,
                               double*
                               >
                             >
                    >
           >::iterator miter1, mend1=sorter.end();
  for (miter1=sorter.begin(); miter1!=mend1; ++miter1)
    {
      std::map<unsigned,
               std::map<unsigned,
                        std::pair<
                          std::pair<double*, double*>,
                          double*
                          >
                        >
               >::iterator miter2, mend2=miter1->second.end();
      for (miter2=miter1->second.begin(); miter2!=mend2; ++miter2)
        {
          std::map<unsigned,
                   std::pair<
                     std::pair<double*, double*>,
                     double*
                     >
                   >::iterator miter3, mend3=miter2->second.end();
          for (miter3=miter2->second.begin(); miter3!=mend3; ++miter3)
            {
              GABAcurrent.push_back(miter3->second.first.first);
              Ca.push_back(miter3->second.first.second);
              eCB.push_back(miter3->second.second);
            }
        }
    }

  // Create the output files...
  std::ostringstream sysCall;
  sysCall<<"mkdir -p "<<directory.c_str()<<";";
  try
    {
      int systemRet = system(sysCall.str().c_str());
      if (systemRet == -1)
        throw;
    }
  catch(...) { };

  std::ostringstream os_GABA, os_Ca, os_eCBproduction, os_eCB;

  int Xdim = (int) mxslice+1;
  int Ydim = (int) mxcol+1;
  int Zdim = (int) mxrow+1;

  if (op_saveGABA)
    {
      os_GABA<<directory<<filePrep<<"DendriteGABA"<<fileApp<<fileExt;
      GABA_file=new std::ofstream(os_GABA.str().c_str(),
                                  std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
      GABA_file->write(reinterpret_cast<char *>(&Xdim), sizeof(Xdim));
      GABA_file->write(reinterpret_cast<char *>(&Ydim), sizeof(Ydim));
      GABA_file->write(reinterpret_cast<char *>(&Zdim), sizeof(Zdim));
    }

  if (op_saveCa)
    {
      os_Ca<<directory<<filePrep<<"DendriteCa"<<fileApp<<fileExt;
      Ca_file=new std::ofstream(os_Ca.str().c_str(),
                                std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
      Ca_file->write(reinterpret_cast<char *>(&Xdim), sizeof(Xdim));
      Ca_file->write(reinterpret_cast<char *>(&Ydim), sizeof(Ydim));
      Ca_file->write(reinterpret_cast<char *>(&Zdim), sizeof(Zdim));
    }

  if (op_saveeCB)
    {
      // Save the eCB production function first
      os_eCBproduction<<directory<<filePrep<<"DendriteeCBproduction"<<fileApp<<fileExt;
      eCBproduction_file=new std::ofstream(os_eCBproduction.str().c_str(),
                                           std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
      double Ca = 0.0;
      float ecb = 0.0;
      for (int i=0; i <= 2000; i++)
        {
          ecb = (float) eCBproduction(Ca);
          eCBproduction_file->write(reinterpret_cast<char *>(&ecb), sizeof(ecb));
          Ca += 1.0 / 1000.0;
        }
      eCBproduction_file->close();
      delete eCBproduction_file;

      // Now the actual eCB file
      os_eCB<<directory<<filePrep<<"DendriteeCB"<<fileApp<<fileExt;
      eCB_file=new std::ofstream(os_eCB.str().c_str(),
                                 std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
      eCB_file->write(reinterpret_cast<char *>(&Xdim), sizeof(Xdim));
      eCB_file->write(reinterpret_cast<char *>(&Ydim), sizeof(Ydim));
      eCB_file->write(reinterpret_cast<char *>(&Zdim), sizeof(Zdim));
    }
}

void DendriteIAFUnitDataCollector::finalize(RNG& rng)
{
  // Close the output files...
  if (op_saveGABA)
    {
      GABA_file->close();
      delete GABA_file;
    }
  if (op_saveCa)
    {
      Ca_file->close();
      delete Ca_file;
    }
  if (op_saveeCB)
    {
      eCB_file->close();
      delete eCB_file;
    }
}

void DendriteIAFUnitDataCollector::dataCollection(Trigger* trigger, NDPairList* ndPairList)
{
  if (op_saveGABA)
    {
      ShallowArray<double*>::iterator iter, end;
      float temp = 0.;
      iter=GABAcurrent.begin();
      end=GABAcurrent.end();
      for (int n=0; iter!=end; ++iter)
        {
          temp = (float) **iter;
          GABA_file->write(reinterpret_cast<char *>(&temp), sizeof(temp));
        }
    }

  if (op_saveCa)
    {
      ShallowArray<double*>::iterator iter, end;
      float temp = 0.;
      iter=Ca.begin();
      end=Ca.end();
      for (int n=0; iter!=end; ++iter)
        {
          temp = (float) **iter;
          Ca_file->write(reinterpret_cast<char *>(&temp), sizeof(temp));
        }
    }

  if (op_saveeCB)
    {
      ShallowArray<double*>::iterator iter, end;
      float temp = 0.;
      iter=eCB.begin();
      end=eCB.end();

      for (int n=0; iter!=end; ++iter)
        {
          temp = (float) **iter;
          eCB_file->write(reinterpret_cast<char *>(&temp), sizeof(temp));
        }
    }
}

void DendriteIAFUnitDataCollector::getNodeIndices(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_DendriteIAFUnitDataCollectorInAttrPSet* CG_inAttrPset, CG_DendriteIAFUnitDataCollectorOutAttrPSet* CG_outAttrPset)
{
  ShallowArray<unsigned,3,2> coords;
  CG_node->getNode()->getNodeCoords(coords);
  assert(coords.size()==3);
  rows.push_back(coords[0]);
  cols.push_back(coords[1]);
  slices.push_back(coords[2]);
}

DendriteIAFUnitDataCollector::DendriteIAFUnitDataCollector()
  : CG_DendriteIAFUnitDataCollector()
{
}

DendriteIAFUnitDataCollector::~DendriteIAFUnitDataCollector()
{
}

void DendriteIAFUnitDataCollector::duplicate(std::auto_ptr<DendriteIAFUnitDataCollector>& dup) const
{
  dup.reset(new DendriteIAFUnitDataCollector(*this));
}

void DendriteIAFUnitDataCollector::duplicate(std::auto_ptr<Variable>& dup) const
{
  dup.reset(new DendriteIAFUnitDataCollector(*this));
}

void DendriteIAFUnitDataCollector::duplicate(std::auto_ptr<CG_DendriteIAFUnitDataCollector>& dup) const
{
  dup.reset(new DendriteIAFUnitDataCollector(*this));
}


double DendriteIAFUnitDataCollector::eCBsigmoid(double Ca)
{
  return 1.0 / ( 1.0 + exp(-eCBprodC * (Ca - eCBprodD)) );
}

double DendriteIAFUnitDataCollector::eCBproduction(double Ca)
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
