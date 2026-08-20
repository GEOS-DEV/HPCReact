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

template< typename REAL_TYPE >
class StoichiometricIonicStrength
{
public:

  template< typename ARRAY_1D_TO_CONST >
  static inline HPCREACT_HOST_DEVICE
  REAL_TYPE
  calculate( ARRAY_1D_TO_CONST const & speciesConcentration,
             ARRAY_1D_TO_CONST const & speciesCharge,
             int const numSpecies )
  {
    return 0.0;
  }

};


} // namespace hpcReact
