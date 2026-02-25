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
 * @brief Convert concentrations to activities and (optionally) map derivatives.
 * @details
 *   When `LOGE_CONCENTRATION == false`, this is a direct pass-through to the activity
 *   model (`a(c)` and `da/dc`).
 *   When `LOGE_CONCENTRATION == true`, the input is interpreted as `log(c)` and the
 *   output as `log(a)`, with derivatives returned as `d log(a) / d log(c)`.
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
  HPCREACT_UNUSED_VAR( activityParams );
  HPCREACT_UNUSED_VAR( speciesConcentration );
  HPCREACT_UNUSED_VAR( activities );
  HPCREACT_UNUSED_VAR( dActivities_dConcentration );
  HPCREACT_UNUSED_VAR( sizeof( INT_TYPE ) );

  static constexpr INDEX_TYPE numSpecies = ACTIVITY_MODEL::Params::numSpecies();

  if constexpr( LOGE_CONCENTRATION )
  {
    REAL_TYPE linearConcentration[numSpecies] = { 0.0 };
    REAL_TYPE linearActivities[numSpecies] = { 0.0 };
    REAL_TYPE dLinearActivities_dLinearConcentration[numSpecies][numSpecies] = {{ 0.0 }};

    for( INDEX_TYPE i = 0; i < numSpecies; ++i )
    {
      linearConcentration[i] = exp( speciesConcentration[i] );
    }

    ACTIVITY_MODEL::calculateActivities( activityParams,
                                         linearConcentration,
                                         linearActivities,
                                         dLinearActivities_dLinearConcentration );

    for( INDEX_TYPE i = 0; i < numSpecies; ++i )
    {
      REAL_TYPE const ai = linearActivities[i];
      activities[i] = log( ai );

      for( INDEX_TYPE j = 0; j < numSpecies; ++j )
      {
        dActivities_dConcentration[i][j] = 0.0;
        if( ai > 1.0e-300 )
        {
          dActivities_dConcentration[i][j] =
            dLinearActivities_dLinearConcentration[i][j] * linearConcentration[j] / ai;
        }
      }
    }
  }
  else
  {
    ACTIVITY_MODEL::calculateActivities( activityParams,
                                         speciesConcentration,
                                         activities,
                                         dActivities_dConcentration );
  }
}

} // namespace hpcReact
