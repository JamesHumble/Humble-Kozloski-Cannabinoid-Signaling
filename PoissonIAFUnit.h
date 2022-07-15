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

#ifndef PoissonIAFUnit_H
#define PoissonIAFUnit_H

#include "Lens.h"
#include "CG_PoissonIAFUnit.h"
#include "rndm.h"

class PoissonIAFUnit : public CG_PoissonIAFUnit
{
 public:
  void update(RNG& rng);
  virtual ~PoissonIAFUnit();
};

#endif
