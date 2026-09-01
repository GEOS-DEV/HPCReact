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

/**
 * @brief The molal ionic strength I = 0.5 * sum_i c_i z_i^2, formed from the speciated
 *        concentrations.
 * @tparam REAL_TYPE floating point type.
 * @tparam INDEX_TYPE integral type used to index the species.
 * @tparam NUM_SPECIES number of species in the system.
 */
template< typename REAL_TYPE,
          typename INDEX_TYPE,
          int NUM_SPECIES >
class SpeciatedIonicStrength
{
public:
  /// alias for the floating point type used in the class.
  using RealType = REAL_TYPE;

  /// alias for the integral type used to index the species.
  using IndexType = INDEX_TYPE;

  /// The parameters the ionic strength requires, and the base of every activity model's Params.
  struct Params
  {

    /// @return The number of species in the system.
    HPCREACT_HOST_DEVICE static constexpr IndexType numSpecies() { return NUM_SPECIES; }

    /// @return A mutable reference to the array of species charges.
    HPCREACT_HOST_DEVICE constexpr CArrayWrapper< RealType, NUM_SPECIES > & speciesCharge() { return m_speciesCharge; }

    /// Charge z_i of each species.
    CArrayWrapper< RealType, NUM_SPECIES > m_speciesCharge;
  };


  /**
   * @brief Compute the ionic strength, and its derivatives wrt the species concentrations.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @tparam ARRAY_1D The type of the array of ionic strength derivatives.
   * @param params the ionic strength parameters
   * @param speciesConcentration linear concentrations c_i
   * @param dIonicStrength_dConcentration [out] dI / d c_i
   * @return the molal ionic strength I
   */
  template< typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D >
  static inline HPCREACT_HOST_DEVICE
  REAL_TYPE
  calculate( Params const & params,
             ARRAY_1D_TO_CONST const & speciesConcentration,
             ARRAY_1D & dIonicStrength_dConcentration )
  {
    REAL_TYPE I = 0.0;
    auto const & numSpecies = params.numSpecies();
    auto const & speciesCharge = params.m_speciesCharge;
    for( int i=0; i<numSpecies; ++i )
    {
      dIonicStrength_dConcentration[i] = 0.5 * speciesCharge[i] * speciesCharge[i];
      I += speciesConcentration[i] * dIonicStrength_dConcentration[i];
    }
    return I;
  }


};


} // namespace hpcReact
