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

/**
 * @brief The stoichiometric ionic strength, formed from the total (unspeciated) concentrations.
 * @tparam REAL_TYPE floating point type.
 * @note Not implemented: calculate() currently returns 0.
 */
template< typename REAL_TYPE >
class StoichiometricIonicStrength
{
public:

  /**
   * @brief Compute the stoichiometric ionic strength.
   * @tparam ARRAY_1D_TO_CONST The type of the arrays of concentrations and charges.
   * @param speciesConcentration linear concentrations c_i
   * @param speciesCharge charge z_i of each species
   * @param numSpecies number of species in the system
   * @return the stoichiometric ionic strength I
   */
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
