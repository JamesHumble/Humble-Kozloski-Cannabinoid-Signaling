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
#include "CleftAstrocyteIAFUnitDataCollector.h"
#include "CG_CleftAstrocyteIAFUnitDataCollector.h"
#include "NodeDescriptor.h"
#include "Node.h"
#include <map>
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <utility>

void CleftAstrocyteIAFUnitDataCollector::initialize(RNG& rng) 
{  
  // Sort pointers by indices, row major
  std::map<unsigned, 
	   std::map<unsigned, 
                    std::map<unsigned,
                             std::pair<
                               std::pair<double*, double*>, // second.first.first second.first.second
                               std::pair<double*, double*>  // second.second.first second.second.second
                               >
                             >
                    >
           >
    sorter;
  assert(rows.size()==slices.size());
  assert(cols.size()==slices.size());
  assert(slices.size()==glutamate.size());
  assert(slices.size()==GABA.size());
  assert(slices.size()==glutamateeCB.size());
  assert(slices.size()==GABAeCB.size());
  int sz=glutamate.size();
  int mxrow=0;
  int mxcol=0;
  for (int j=0; j<sz; ++j)
    {
      sorter[rows[j]][cols[j]][slices[j]]=std::make_pair(
                                                         std::make_pair(glutamate[j], GABA[j]),
                                                         std::make_pair(glutamateeCB[j], GABAeCB[j])
                                                         );
      if (mxrow<rows[j]) mxrow=rows[j];
      if (mxcol<cols[j]) mxcol=cols[j];
      if (mxslice<slices[j]) mxslice=slices[j];
    }
  glutamate.clear();
  GABA.clear();
  glutamateeCB.clear();
  GABAeCB.clear();
  std::map<unsigned, 
	   std::map<unsigned, 
                    std::map<unsigned,
                             std::pair<
                               std::pair<double*, double*>,
                               std::pair<double*, double*>
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
                          std::pair<double*, double*>
                          >
                        >
               >::iterator miter2, mend2=miter1->second.end();    
      for (miter2=miter1->second.begin(); miter2!=mend2; ++miter2)
        {
          std::map<unsigned,
                   std::pair<
                     std::pair<double*, double*>,
                     std::pair<double*, double*>
                     >
                   >::iterator miter3, mend3=miter2->second.end();
          for (miter3=miter2->second.begin(); miter3!=mend3; ++miter3)
            {
              glutamate.push_back(miter3->second.first.first);
              GABA.push_back(miter3->second.first.second);
              glutamateeCB.push_back(miter3->second.second.first);
              GABAeCB.push_back(miter3->second.second.second);
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
  
  std::ostringstream os_glutamate, os_GABA, os_glutamateeCB, os_GABAeCB;

  int Xdim = (int) mxslice+1;
  int Ydim = (int) mxcol+1;
  int Zdim = (int) mxrow+1;

  if (op_saveGlutamate)
    {
      os_glutamate<<directory<<filePrep<<"CleftAstrocyteGlutamate"
                         <<fileApp<<fileExt;
      glutamate_file=new std::ofstream(os_glutamate.str().c_str(),
                                       std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
      glutamate_file->write(reinterpret_cast<char *>(&Xdim), sizeof(Xdim));
      glutamate_file->write(reinterpret_cast<char *>(&Ydim), sizeof(Ydim));
      glutamate_file->write(reinterpret_cast<char *>(&Zdim), sizeof(Zdim));
    }

  if (op_saveGABA)
    {
      os_GABA<<directory<<filePrep<<"CleftAstrocyteGABA"
                         <<fileApp<<fileExt;
      GABA_file=new std::ofstream(os_GABA.str().c_str(),
                                       std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
      GABA_file->write(reinterpret_cast<char *>(&Xdim), sizeof(Xdim));
      GABA_file->write(reinterpret_cast<char *>(&Ydim), sizeof(Ydim));
      GABA_file->write(reinterpret_cast<char *>(&Zdim), sizeof(Zdim));
    }  

  if (op_saveGlutamateeCB)
    {
      os_glutamateeCB<<directory<<filePrep<<"CleftAstrocyteGlutamateeCB"<<fileApp<<fileExt;
      glutamateeCB_file=new std::ofstream(os_glutamateeCB.str().c_str(),
                                 std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
      glutamateeCB_file->write(reinterpret_cast<char *>(&Xdim), sizeof(Xdim));
      glutamateeCB_file->write(reinterpret_cast<char *>(&Ydim), sizeof(Ydim));
      glutamateeCB_file->write(reinterpret_cast<char *>(&Zdim), sizeof(Zdim));
    }

  if (op_saveGABAeCB)
    {
      os_GABAeCB<<directory<<filePrep<<"CleftAstrocyteGABAeCB"<<fileApp<<fileExt;
      GABAeCB_file=new std::ofstream(os_GABAeCB.str().c_str(),
                                 std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
      GABAeCB_file->write(reinterpret_cast<char *>(&Xdim), sizeof(Xdim));
      GABAeCB_file->write(reinterpret_cast<char *>(&Ydim), sizeof(Ydim));
      GABAeCB_file->write(reinterpret_cast<char *>(&Zdim), sizeof(Zdim));
    }        
}

void CleftAstrocyteIAFUnitDataCollector::finalize(RNG& rng) 
{
  // Close the output files...
  if (op_saveGlutamate)
    {
      glutamate_file->close();
      delete glutamate_file;
    }
  if (op_saveGABA)
    {
      GABA_file->close();
      delete GABA_file;
    }    
  if (op_saveGlutamateeCB)
    {
      glutamateeCB_file->close();
      delete glutamateeCB_file;
    }
  if (op_saveGABAeCB)
    {
      GABAeCB_file->close();
      delete GABAeCB_file;
    }    
}

void CleftAstrocyteIAFUnitDataCollector::dataCollection(Trigger* trigger, NDPairList* ndPairList) 
{
  if (op_saveGlutamate)
    {
      ShallowArray<double*>::iterator iter=glutamate.begin(), end=glutamate.end();
      float temp = 0.;
      for (int n=0; iter!=end; ++iter, n++)
        {
          temp = (float) **iter;
          glutamate_file->write(reinterpret_cast<char *>(&temp), sizeof(temp));            
        }
    }

  if (op_saveGABA)
    {
      ShallowArray<double*>::iterator iter=GABA.begin(), end=GABA.end();
      float temp = 0.;
      for (int n=0; iter!=end; ++iter, n++)
        {
          temp = (float) **iter;
          GABA_file->write(reinterpret_cast<char *>(&temp), sizeof(temp));            
        }
    }  

  if (op_saveGlutamateeCB)
    {
      ShallowArray<double*>::iterator iter=glutamateeCB.begin(), end=glutamateeCB.end();
      float temp = 0.;
      for (int n=0; iter!=end; ++iter, n++)
        {
          temp = (float) **iter;
          glutamateeCB_file->write(reinterpret_cast<char *>(&temp), sizeof(temp));            
        }
    }

  if (op_saveGABAeCB)
    {
      ShallowArray<double*>::iterator iter=GABAeCB.begin(), end=GABAeCB.end();
      float temp = 0.;
      for (int n=0; iter!=end; ++iter, n++)
        {
          temp = (float) **iter;
          GABAeCB_file->write(reinterpret_cast<char *>(&temp), sizeof(temp));            
        }
    }    
}

void CleftAstrocyteIAFUnitDataCollector::getNodeIndices(const String& CG_direction, const String& CG_component, NodeDescriptor* CG_node, Edge* CG_edge, VariableDescriptor* CG_variable, Constant* CG_constant, CG_CleftAstrocyteIAFUnitDataCollectorInAttrPSet* CG_inAttrPset, CG_CleftAstrocyteIAFUnitDataCollectorOutAttrPSet* CG_outAttrPset) 
{
  ShallowArray<unsigned,3,2> coords;
  CG_node->getNode()->getNodeCoords(coords);
  assert(coords.size()==3);
  rows.push_back(coords[0]);
  cols.push_back(coords[1]);
  slices.push_back(coords[2]);
}

CleftAstrocyteIAFUnitDataCollector::CleftAstrocyteIAFUnitDataCollector() 
  : CG_CleftAstrocyteIAFUnitDataCollector()
{
}

CleftAstrocyteIAFUnitDataCollector::~CleftAstrocyteIAFUnitDataCollector() 
{
}

void CleftAstrocyteIAFUnitDataCollector::duplicate(std::auto_ptr<CleftAstrocyteIAFUnitDataCollector>& dup) const
{
  dup.reset(new CleftAstrocyteIAFUnitDataCollector(*this));
}

void CleftAstrocyteIAFUnitDataCollector::duplicate(std::auto_ptr<Variable>& dup) const
{
  dup.reset(new CleftAstrocyteIAFUnitDataCollector(*this));
}

void CleftAstrocyteIAFUnitDataCollector::duplicate(std::auto_ptr<CG_CleftAstrocyteIAFUnitDataCollector>& dup) const
{
  dup.reset(new CleftAstrocyteIAFUnitDataCollector(*this));
}
