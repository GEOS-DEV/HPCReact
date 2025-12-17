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
#include "reactions/reactionsSystems/MixedEquilibriumKineticReactions.hpp"
#include "reactions/reactionsSystems/EquilibriumReactions.hpp"
#include "common/macros.hpp"
#include "common/printers.hpp"
#include "common/nonlinearSolvers.hpp"
#include "common/pmpl.hpp"

#include <gtest/gtest.h>

namespace hpcReact
{

namespace unitTest_utilities
{

//******************************************************************************
template< typename REAL_TYPE,
          typename PARAMS_DATA >
void timeStepTest( PARAMS_DATA const & params,
                   REAL_TYPE const dt,
                   int const numSteps,
                   REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numPrimarySpecies()],
                   REAL_TYPE const (&surfaceArea)[PARAMS_DATA::numKineticReactions()],
                   REAL_TYPE const (&expectedSpeciesConcentrations)[PARAMS_DATA::numPrimarySpecies()] )
{
  HPCREACT_UNUSED_VAR( expectedSpeciesConcentrations );

  CArrayWrapper< REAL_TYPE, PARAMS_DATA::numPrimarySpecies() > primarySpeciesConcentration;

  for( int i = 0; i < PARAMS_DATA::numPrimarySpecies(); ++i )
  {
    primarySpeciesConcentration[i] = initialSpeciesConcentration[i];
  }

  pmpl::genericKernelWrapper( PARAMS_DATA::numPrimarySpecies(),
                              primarySpeciesConcentration.data,
                              [=] HPCREACT_HOST_DEVICE ( decltype(primarySpeciesConcentration.data) speciesConcentration )
      {
        using MixedReactionsType = reactionsSystems::MixedEquilibriumKineticReactions< REAL_TYPE,
                                                                                       int,
                                                                                       int >;
        using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< REAL_TYPE,
                                                                                 int,
                                                                                 int >;

        // constexpr int numSpecies = PARAMS_DATA::numSpecies();
        static constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();
        static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
        static constexpr int numKineticReactions = PARAMS_DATA::numKineticReactions();

        // define variables
        double const temperature = 298.15;
        REAL_TYPE logPrimarySpeciesConcentration[numPrimarySpecies];
        // must use CArrayWrapper to ensure correct capture in the lambda functions
        CArrayWrapper< REAL_TYPE, numSecondarySpecies > logSecondarySpeciesConcentration;
        CArrayWrapper< REAL_TYPE, numPrimarySpecies > aggregatePrimarySpeciesConcentration;
        CArrayWrapper< REAL_TYPE, numPrimarySpecies > aggregatePrimarySpeciesConcentration_n;
        CArrayWrapper< REAL_TYPE, numPrimarySpecies > mobileAggregatePrimarySpeciesConcentration;
        CArrayWrapper< REAL_TYPE, numPrimarySpecies, numPrimarySpecies > dAggregatePrimarySpeciesConcentrations_dlogPrimarySpeciesConcentration;
        CArrayWrapper< REAL_TYPE, numPrimarySpecies, numPrimarySpecies > dMobileAggregatePrimarySpeciesConcentrations_dlogPrimarySpeciesConcentration;
        CArrayWrapper< REAL_TYPE, numKineticReactions > reactionRates;
        CArrayWrapper< REAL_TYPE, numKineticReactions, numPrimarySpecies > dReactionRates_dlogPrimarySpeciesConcentration;
        CArrayWrapper< REAL_TYPE, numPrimarySpecies > aggregateSpeciesRates;
        CArrayWrapper< REAL_TYPE, numPrimarySpecies, numPrimarySpecies > dAggregateSpeciesRates_dlogPrimarySpeciesConcentration;

        // Initialize species concentrations
        for( int i = 0; i < numPrimarySpecies; ++i )
        {
          logPrimarySpeciesConcentration[i] = logmath::log( speciesConcentration[i] );
          aggregatePrimarySpeciesConcentration[i] = speciesConcentration[i];
        }

        EquilibriumReactionsType::enforceEquilibrium_LogAggregate( temperature,
                                                                   params.equilibriumReactionsParameters(),
                                                                   logPrimarySpeciesConcentration,
                                                                   logPrimarySpeciesConcentration );

        /// Time step loop
        double time = 0.0;
        for( int t = 0; t < numSteps; ++t )
        {
          printf( "Timestep %d, Time = %.6e\n", t, time );

          for( int i=0; i < numPrimarySpecies; ++i )
          {
            aggregatePrimarySpeciesConcentration_n[i] = aggregatePrimarySpeciesConcentration[i];
          }

          auto computeResidualAndJacobian = [&] ( REAL_TYPE const (&X)[numPrimarySpecies],
                                                                       REAL_TYPE ( &r )[numPrimarySpecies],
                                                                       REAL_TYPE ( &J )[numPrimarySpecies][numPrimarySpecies] )
          {
            MixedReactionsType::updateMixedSystem( temperature,
                                                   params,
                                                   X,
                                                   surfaceArea,
                                                   logSecondarySpeciesConcentration,
                                                   aggregatePrimarySpeciesConcentration,
                                                   mobileAggregatePrimarySpeciesConcentration,
                                                   dAggregatePrimarySpeciesConcentrations_dlogPrimarySpeciesConcentration,
                                                   dMobileAggregatePrimarySpeciesConcentrations_dlogPrimarySpeciesConcentration,
                                                   reactionRates,
                                                   dReactionRates_dlogPrimarySpeciesConcentration,
                                                   aggregateSpeciesRates,
                                                   dAggregateSpeciesRates_dlogPrimarySpeciesConcentration );


            for( int i = 0; i < numPrimarySpecies; ++i )
            {
              r[i] = ( aggregatePrimarySpeciesConcentration[i] - aggregatePrimarySpeciesConcentration_n[i] ) - aggregateSpeciesRates[i] * dt;
              for( int j = 0; j < numPrimarySpecies; ++j )
              {
                J[i][j] = dAggregatePrimarySpeciesConcentrations_dlogPrimarySpeciesConcentration[i][j] - dAggregateSpeciesRates_dlogPrimarySpeciesConcentration[i][j] * dt;
              }
            }
          };

          nonlinearSolvers::newtonRaphson< numPrimarySpecies >( logPrimarySpeciesConcentration, computeResidualAndJacobian );

          time += dt;
        }
        for( int i = 0; i < numPrimarySpecies; ++i )
        {
          speciesConcentration[i] = logmath::exp( logPrimarySpeciesConcentration[i] );
        }
      } );

  // Check results
  for( int i = 0; i < PARAMS_DATA::numPrimarySpecies(); ++i )
  {
    EXPECT_NEAR( primarySpeciesConcentration[ i ], expectedSpeciesConcentrations[ i ], 1.0e-8 * expectedSpeciesConcentrations[ i ] );
  }
}


} // namespace unitTest_utilities
} // namespace hpcReact
