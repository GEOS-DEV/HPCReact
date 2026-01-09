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
class Identity
{
public:
  using RealType = REAL_TYPE;
  using IndexType = INDEX_TYPE;

  struct Params
  { 
    HPCREACT_HOST_DEVICE static constexpr IndexType numSpecies() { return IONIC_STRENGTH_TYPE::Params::numSpecies(); }
  };

  template< typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename PARAMS >
  static inline HPCREACT_HOST_DEVICE
  void
  calculateActivities( PARAMS const & ,
                       ARRAY_1D_TO_CONST const & speciesConcentrations,
                       ARRAY_1D & activities )
  {
    constexpr IndexType numSpecies = PARAMS::numSpecies();
    for( IndexType i=0; i<numSpecies; ++i )
    {
      activities[i] = speciesConcentrations[i];
    }
  }

};


} // namespace hpcReact
