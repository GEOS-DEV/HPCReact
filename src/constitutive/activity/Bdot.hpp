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

#include "common/macros.hpp"

namespace hpcReact
{

template< typename REAL_TYPE,
          typename INDEX_TYPE,
          typename IONIC_STRENGTH_TYPE >
class Bdot
{
public:
  using RealType = REAL_TYPE;
  using IndexType = INDEX_TYPE;

struct Params : public IONIC_STRENGTH_TYPE::PARAMS
{

};


template< typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D >
static inline HPCREACT_HOST_DEVICE 
void
calculateActivities( IONIC_STRENGTH_TYPE::PARAMS const & params,
                     ARRAY_1D_TO_CONST const & speciesConcentrations,
                     ARRAY_1D & activities )
{

  RealType const ionicStrength = IONIC_STRENGTH_TYPE::calculate( params, speciesConcentrations );
  RealType const sqrtI = sqrt(ionicStrength);
  RealType const A_gamma = 2;
  RealType const B_gamma = 1.6;
  auto const & speciesCharge = params.m_speciesCharge;
  auto const & a = params.m_ionSizeParameter;
  auto const & b = params.m_bdotParameter;



  constexpr IndexType numSpecies = params.numSpecies();
  for( IndexType i=0; i<numSpecies; ++i )
  {
    RealType const gamma_coeff = -A_gamma * sqrtI / ( 1 + a[i] * B_gamma * sqrtI );
    activities[i] = gamma_coeff * speciesCharge[i] * speciesCharge[i] + b[i] * ionicStrength;
  }
}

};


} // namespace hpcReact