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

#include <math.h>

namespace hpcReact
{

/**
 * @brief Convert concentrations to activities, exposing the activity coefficients used.
 * @details
 *   The model supplies `ln(gamma)`; this composes `a = c * gamma`.
 *   `LOGE_CONCENTRATION == false`: input `c`, output `a`, derivatives wrt `c`.
 *   `LOGE_CONCENTRATION == true`: input `log(c)`, output `log(a)`, derivatives wrt `log(c)`.
 *   `logActivityCoefficients` is `ln(gamma)` in either case, for callers that need `C = a / gamma`.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL,
          bool LOGE_CONCENTRATION,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D,
          typename ARRAY_1D_GAMMA,
          typename ARRAY_2D_GAMMA >
HPCREACT_HOST_DEVICE
inline
void calculateActivities( typename ACTIVITY_MODEL::Params const & activityParams,
                          ARRAY_1D_TO_CONST const & speciesConcentration,
                          ARRAY_1D & activities,
                          ARRAY_2D & dActivities_dConcentration,
                          ARRAY_1D_GAMMA & logActivityCoefficients,
                          ARRAY_2D_GAMMA & dLogActivityCoefficients_dConcentration )
{
  HPCREACT_UNUSED_VAR( sizeof( INT_TYPE ) );

  static constexpr INDEX_TYPE numSpecies = ACTIVITY_MODEL::Params::numSpecies();

  if constexpr( LOGE_CONCENTRATION )
  {
    REAL_TYPE linearConcentration[numSpecies] = { 0.0 };

    for( INDEX_TYPE i = 0; i < numSpecies; ++i )
    {
      linearConcentration[i] = exp( speciesConcentration[i] );
    }

    ACTIVITY_MODEL::calculateLogActivityCoefficients( activityParams,
                                                      linearConcentration,
                                                      logActivityCoefficients,
                                                      dLogActivityCoefficients_dConcentration );

    for( INDEX_TYPE i = 0; i < numSpecies; ++i )
    {
      activities[i] = speciesConcentration[i] + logActivityCoefficients[i];

      for( INDEX_TYPE j = 0; j < numSpecies; ++j )
      {
        // d ln(gamma_i)/d ln(c_j) = d ln(gamma_i)/d c_j * c_j
        REAL_TYPE const dLogGamma_dLogC =
          dLogActivityCoefficients_dConcentration[i][j] * linearConcentration[j];

        dLogActivityCoefficients_dConcentration[i][j] = dLogGamma_dLogC;
        dActivities_dConcentration[i][j] = ( i == j ? 1.0 : 0.0 ) + dLogGamma_dLogC;
      }
    }
  }
  else
  {
    ACTIVITY_MODEL::calculateLogActivityCoefficients( activityParams,
                                                      speciesConcentration,
                                                      logActivityCoefficients,
                                                      dLogActivityCoefficients_dConcentration );

    for( INDEX_TYPE i = 0; i < numSpecies; ++i )
    {
      REAL_TYPE const gamma_i = exp( logActivityCoefficients[i] );
      activities[i] = speciesConcentration[i] * gamma_i;

      for( INDEX_TYPE j = 0; j < numSpecies; ++j )
      {
        dActivities_dConcentration[i][j] =
          ( i == j ? gamma_i : 0.0 )
          + speciesConcentration[i] * gamma_i * dLogActivityCoefficients_dConcentration[i][j];
      }
    }
  }
}

/**
 * @brief Convert concentrations to activities, discarding the activity coefficients.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL,
          bool LOGE_CONCENTRATION,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateActivities( typename ACTIVITY_MODEL::Params const & activityParams,
                          ARRAY_1D_TO_CONST const & speciesConcentration,
                          ARRAY_1D & activities,
                          ARRAY_2D & dActivities_dConcentration )
{
  static constexpr INDEX_TYPE numSpecies = ACTIVITY_MODEL::Params::numSpecies();

  REAL_TYPE logActivityCoefficients[numSpecies] = { 0.0 };
  REAL_TYPE dLogActivityCoefficients_dConcentration[numSpecies][numSpecies] = {{ 0.0 }};

  calculateActivities< REAL_TYPE,
                       INT_TYPE,
                       INDEX_TYPE,
                       ACTIVITY_MODEL,
                       LOGE_CONCENTRATION >( activityParams,
                                             speciesConcentration,
                                             activities,
                                             dActivities_dConcentration,
                                             logActivityCoefficients,
                                             dLogActivityCoefficients_dConcentration );
}

} // namespace hpcReact
