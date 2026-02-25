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

#include "common/constants.hpp"
#include "common/CArrayWrapper.hpp"
#include "reactions/massActions/MassActions.hpp"


/** @file MixedEquilibriumKineticReactions_impl.hpp
 *  @brief Header file for the MixedEquilibriumKineticReactions implementation.
 *  @author HPC-REACT Team
 *  @date 2025
 */

namespace hpcReact
{
namespace reactionsSystems
{

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL,
          bool LOGE_CONCENTRATION >
template< typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST_KINETIC,
          typename ARRAY_1D_PRIMARY,
          typename ARRAY_1D_SECONDARY,
          typename ARRAY_1D_KINETIC,
          typename ARRAY_2D_PRIMARY,
          typename ARRAY_2D_KINETIC >
HPCREACT_HOST_DEVICE inline void
MixedEquilibriumKineticReactions< REAL_TYPE,
                                  INT_TYPE,
                                  INDEX_TYPE,
                                  ACTIVITY_MODEL,
                                  LOGE_CONCENTRATION
                                  >::updateMixedSystem_impl( RealType const & temperature,
                                                             PARAMS_DATA const & params,
                                                             typename ACTIVITY_MODEL::Params const & activityParams,
                                                             ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                             ARRAY_1D_TO_CONST_KINETIC const & surfaceArea,
                                                             ARRAY_1D_SECONDARY & logSecondarySpeciesConcentrations,
                                                             ARRAY_1D_PRIMARY & aggregatePrimarySpeciesConcentrations,
                                                             ARRAY_1D_PRIMARY & mobileAggregatePrimarySpeciesConcentrations,
                                                             ARRAY_2D_PRIMARY & dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                                                             ARRAY_2D_PRIMARY & dMobileAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                                                             ARRAY_1D_KINETIC & reactionRates,
                                                             ARRAY_2D_KINETIC & dReactionRates_dLogPrimarySpeciesConcentrations,
                                                             ARRAY_1D_PRIMARY & aggregateSpeciesRates,
                                                             ARRAY_2D_PRIMARY & dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations )
{
  constexpr IntType numSpecies = PARAMS_DATA::numSpecies();
  constexpr IntType numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  constexpr IntType numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  if constexpr( PARAMS_DATA::numEquilibriumReactions() > 0 )
  {
    RealType logSpeciesConcentration[numSpecies] = { 0.0 };
    RealType logSpeciesActivities[numSpecies] = { 0.0 };
    RealType dLogSpeciesActivities_dLogSpeciesConcentrations[numSpecies][numSpecies] = {{ 0.0 }};
    RealType logPrimaryActivities[numPrimarySpecies] = { 0.0 };
    RealType dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[numPrimarySpecies][numPrimarySpecies] = {{ 0.0 }};
    RealType logSecondarySpeciesConcentrations_guess[numSecondarySpecies] = { 0.0 };
    RealType dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[numSecondarySpecies][numPrimarySpecies] = {{ 0.0 }};

    massActions::calculateLogSecondarySpeciesConcentrationWrtLogC< REAL_TYPE,
                                                                   INT_TYPE,
                                                                   INDEX_TYPE >( params.equilibriumReactionsParameters(),
                                                                                 logPrimarySpeciesConcentrations,
                                                                                 logSecondarySpeciesConcentrations_guess,
                                                                                 dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations );
    for( IntType i = 0; i < numSecondarySpecies; ++i )
    {
      logSpeciesConcentration[i] = logSecondarySpeciesConcentrations_guess[i];
    }
    for( IntType i = 0; i < numPrimarySpecies; ++i )
    {
      logSpeciesConcentration[i + numSecondarySpecies] = logPrimarySpeciesConcentrations[i];
    }

    calculateActivities< RealType,
                         IntType,
                         IndexType,
                         ACTIVITY_MODEL,
                         true >( activityParams,
                                 logSpeciesConcentration,
                                 logSpeciesActivities,
                                 dLogSpeciesActivities_dLogSpeciesConcentrations );

    for( IntType i = 0; i < numPrimarySpecies; ++i )
    {
      IntType const fullRow = i + numSecondarySpecies;
      logPrimaryActivities[i] = logSpeciesActivities[fullRow];

      for( IntType j = 0; j < numPrimarySpecies; ++j )
      {
        RealType value = dLogSpeciesActivities_dLogSpeciesConcentrations[fullRow][j + numSecondarySpecies];
        for( IntType k = 0; k < numSecondarySpecies; ++k )
        {
          value += dLogSpeciesActivities_dLogSpeciesConcentrations[fullRow][k] *
                   dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[k][j];
        }
        dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[i][j] = value;
      }
    }

    // 1. Compute new aggregate species from primary species
    massActions::calculateTotalAndMobileAggregatePrimaryConcentrationsWrtLogC< REAL_TYPE,
                                                                               INT_TYPE,
                                                                               INDEX_TYPE >( params.equilibriumReactionsParameters(),
                                                                                             logPrimarySpeciesConcentrations,
                                                                                             logPrimaryActivities,
                                                                                             dLogPrimaryActivities_dLogPrimarySpeciesConcentrations,
                                                                                             logSecondarySpeciesConcentrations,
                                                                                             aggregatePrimarySpeciesConcentrations,
                                                                                             mobileAggregatePrimarySpeciesConcentrations,
                                                                                             dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                                                                                             dMobileAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations );
  }
  else
  {
    for( int i = 0; i < numPrimarySpecies; ++i )
    {
      REAL_TYPE const speciesConcentration_i = exp( logPrimarySpeciesConcentrations[i] );
      aggregatePrimarySpeciesConcentrations[i] = speciesConcentration_i;
      mobileAggregatePrimarySpeciesConcentrations[i] = speciesConcentration_i;
      dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations( i, i ) = speciesConcentration_i;
      dMobileAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations( i, i ) = speciesConcentration_i;
    }
  }

  if constexpr( PARAMS_DATA::numKineticReactions() > 0 )
  {
    // 2. Compute the reaction rates for all kinetic reactions
    computeReactionRates( temperature,
                          params,
                          activityParams,
                          logPrimarySpeciesConcentrations,
                          logSecondarySpeciesConcentrations,
                          surfaceArea,
                          reactionRates,
                          dReactionRates_dLogPrimarySpeciesConcentrations );

    // 3. Compute aggregate species rates
    computeAggregateSpeciesRates( params,
                                  logPrimarySpeciesConcentrations,
                                  reactionRates,
                                  dReactionRates_dLogPrimarySpeciesConcentrations,
                                  aggregateSpeciesRates,
                                  dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations );
  }
  else
  {
    HPCREACT_UNUSED_VAR( reactionRates );
    HPCREACT_UNUSED_VAR( dReactionRates_dLogPrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( aggregateSpeciesRates );
    HPCREACT_UNUSED_VAR( dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations );
  }

}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL,
          bool LOGE_CONCENTRATION >
template< typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_1D_TO_CONST_KINETIC,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline void
MixedEquilibriumKineticReactions< REAL_TYPE,
                                  INT_TYPE,
                                  INDEX_TYPE,
                                  ACTIVITY_MODEL,
                                  LOGE_CONCENTRATION
                                  >::computeReactionRates_impl( RealType const & temperature,
                                                                PARAMS_DATA const & params,
                                                                typename ACTIVITY_MODEL::Params const & activityParams,
                                                                ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                                ARRAY_1D_TO_CONST2 const & logSecondarySpeciesConcentrations,
                                                                ARRAY_1D_TO_CONST_KINETIC const & surfaceArea,
                                                                ARRAY_1D & reactionRates,
                                                                ARRAY_2D & dReactionRates_dLogPrimarySpeciesConcentrations )

{
  constexpr IntType numSpecies          = PARAMS_DATA::numSpecies();
  constexpr IntType numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  constexpr IntType numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();
  constexpr IntType numKineticReactions = PARAMS_DATA::numKineticReactions();

  RealType logSpeciesConcentration[numSpecies] {};
  for( INDEX_TYPE i = 0; i < numSecondarySpecies; ++i )
  {
    logSpeciesConcentration[i] = logSecondarySpeciesConcentrations[i];
  }
  for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
  {
    logSpeciesConcentration[i+numSecondarySpecies] = logPrimarySpeciesConcentrations[i];
  }

  for( INDEX_TYPE i = 0; i < numKineticReactions; ++i )
  {
    for( INDEX_TYPE j = 0; j < numPrimarySpecies; ++j )
    {
      dReactionRates_dLogPrimarySpeciesConcentrations[i][j] = 0.0;
    }
  }

  CArrayWrapper< RealType, numKineticReactions, numSpecies > reactionRatesDerivatives;

  kineticReactions::computeReactionRates( temperature,
                                          params.kineticReactionsParameters(),
                                          activityParams,
                                          logSpeciesConcentration,
                                          surfaceArea,
                                          reactionRates,
                                          reactionRatesDerivatives );

  // Compute the reaction rates derivatives w.r.t. log primary species concentrations
  for( IntType i = 0; i < numKineticReactions; ++i )
  {
    for( IntType j = 0; j < numPrimarySpecies; ++j )
    {
      dReactionRates_dLogPrimarySpeciesConcentrations( i, j ) = reactionRatesDerivatives( i, j + numSecondarySpecies );
    }
  }

  if constexpr( numSecondarySpecies > 0 )
  {
    RealType logSpeciesActivities[numSpecies] = { 0.0 };
    RealType dLogSpeciesActivities_dLogSpeciesConcentrations[numSpecies][numSpecies] = {{ 0.0 }};
    RealType logPrimaryActivities[numPrimarySpecies] = { 0.0 };
    RealType dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[numPrimarySpecies][numPrimarySpecies] = {{ 0.0 }};
    RealType logSecondaryActivities[numSecondarySpecies] = { 0.0 };
    RealType dLogSecondaryActivities_dLogPrimaryActivities[numSecondarySpecies][numPrimarySpecies] = {{ 0.0 }};
    RealType dLogSecondaryActivities_dLogPrimarySpeciesConcentrations[numSecondarySpecies][numPrimarySpecies] = {{ 0.0 }};

    calculateActivities< RealType,
                         IntType,
                         IndexType,
                         ACTIVITY_MODEL,
                         true >( activityParams,
                                 logSpeciesConcentration,
                                 logSpeciesActivities,
                                 dLogSpeciesActivities_dLogSpeciesConcentrations );

    for( IntType i = 0; i < numPrimarySpecies; ++i )
    {
      IntType const fullRow = i + numSecondarySpecies;
      logPrimaryActivities[i] = logSpeciesActivities[fullRow];
      for( IntType j = 0; j < numPrimarySpecies; ++j )
      {
        dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[i][j] =
          dLogSpeciesActivities_dLogSpeciesConcentrations[fullRow][j + numSecondarySpecies];
      }
    }

    massActions::calculateLogSecondaryActivitiesWrtLogC< REAL_TYPE,
                                                         INT_TYPE,
                                                         INDEX_TYPE >( params.equilibriumReactionsParameters(),
                                                                       logPrimaryActivities,
                                                                       logSecondaryActivities,
                                                                       dLogSecondaryActivities_dLogPrimaryActivities );
    HPCREACT_UNUSED_VAR( logSecondaryActivities );

    for( IntType k = 0; k < numSecondarySpecies; ++k )
    {
      for( IntType j = 0; j < numPrimarySpecies; ++j )
      {
        dLogSecondaryActivities_dLogPrimarySpeciesConcentrations[k][j] = 0.0;
        for( IntType l = 0; l < numPrimarySpecies; ++l )
        {
          dLogSecondaryActivities_dLogPrimarySpeciesConcentrations[k][j] +=
            dLogSecondaryActivities_dLogPrimaryActivities[k][l] *
            dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[l][j];
        }
      }
    }

    for( IntType i = 0; i < numKineticReactions; ++i )
    {
      for( IntType j = 0; j < numPrimarySpecies; ++j )
      {
        for( IntType k = 0; k < numSecondarySpecies; ++k )
        {
          dReactionRates_dLogPrimarySpeciesConcentrations( i, j ) +=
            reactionRatesDerivatives( i, k ) * dLogSecondaryActivities_dLogPrimarySpeciesConcentrations[k][j];
        }
      }
    }
  }
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL,
          bool LOGE_CONCENTRATION >
template< typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D,
          bool CALCULATE_DERIVATIVES >
HPCREACT_HOST_DEVICE inline void
MixedEquilibriumKineticReactions< REAL_TYPE,
                                  INT_TYPE,
                                  INDEX_TYPE,
                                  ACTIVITY_MODEL,
                                  LOGE_CONCENTRATION
                                  >::computeAggregateSpeciesRates_impl( PARAMS_DATA const & params,
                                                                        ARRAY_1D_TO_CONST const & speciesConcentration,
                                                                        ARRAY_1D_TO_CONST2 const & reactionRates,
                                                                        ARRAY_2D_TO_CONST const & reactionRatesDerivatives,
                                                                        ARRAY_1D & aggregatesRates,
                                                                        ARRAY_2D & aggregatesRatesDerivatives )
{
  HPCREACT_UNUSED_VAR( speciesConcentration );
  // constexpr IntType numSpecies = PARAMS_DATA::numSpecies();
  constexpr IntType numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  constexpr IntType numKineticReactions = PARAMS_DATA::numKineticReactions();
  constexpr IntType numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
  {
    aggregatesRates[i] = 0.0;
    if constexpr( CALCULATE_DERIVATIVES )
    {
      for( INDEX_TYPE j = 0; j < numPrimarySpecies; ++j )
      {
        aggregatesRatesDerivatives( i, j ) = 0.0;
      }
    }
    for( IntType r=0; r<numKineticReactions; ++r )
    {
      RealType const s_ir = params.kineticReactionsParameters().stoichiometricMatrix( r, i+numSecondarySpecies );
      aggregatesRates[i] += s_ir * reactionRates[r];
      if constexpr( CALCULATE_DERIVATIVES )
      {
        for( IntType j = 0; j < numPrimarySpecies; ++j )
        {
          aggregatesRatesDerivatives( i, j ) += s_ir * reactionRatesDerivatives( r, j );
        }
      }
    }
  }

}

} // namespace reactionsSystems

} // namespace hpcReact

#include "common/macrosCleanup.hpp"
