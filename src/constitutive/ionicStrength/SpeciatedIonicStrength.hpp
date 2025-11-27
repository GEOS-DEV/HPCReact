/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: (BSD-3-Clause)
 *
 * Copyright (c) 2025- Lawrence Livermore National Security LLC
 * All rights reserved
 *
 * See top level LICENSE files for details.
 * ------------------------------------------------------------------------------------------------------------
 */
#pragma once

#include "common/CArrayWrapper.hpp"
#include "common/macros.hpp"

namespace hpcReact
{

template< typename REAL_TYPE,
          typename INDEX_TYPE,
          int NUM_SPECIES >
class SpeciatedIonicStrength
{
public:
  using RealType = REAL_TYPE;
  using IndexType = INDEX_TYPE;

struct Params
{
  
  HPCREACT_HOST_DEVICE static constexpr IndexType numSpecies() { return NUM_SPECIES; }
  HPCREACT_HOST_DEVICE constexpr CArrayWrapper< RealType, NUM_SPECIES > & speciesCharge() { return m_speciesCharge; }

  CArrayWrapper< RealType, NUM_SPECIES > m_speciesCharge;
};


template< typename ARRAY_1D_TO_CONST >
static inline HPCREACT_HOST_DEVICE 
REAL_TYPE
calculate( Params const & params,
           ARRAY_1D_TO_CONST const & speciesConcentration )
{
  REAL_TYPE I = 0.0;
  auto const &  numSpecies = params.numSpecies();
  auto const & speciesCharge = params.m_speciesCharge;
  for( int i=0; i<numSpecies; ++i )
  {
    I += speciesConcentration[i] * speciesCharge[i] * speciesCharge[i];
  }
  I *= 0.5;
  return I;
}


};


} // namespace hpcReact