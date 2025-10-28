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
#include "common/DirectSystemSolve.hpp"
#include "common/macros.hpp"
#include "common/pmpl.hpp"
#include "common/printers.hpp"
#include "common/CArrayWrapper.hpp"
#include "reactions/reactionsSystems/KineticReactions.hpp"

#include <gtest/gtest.h>

namespace hpcReact
{

namespace unitTest_utilities
{

template< typename REAL_TYPE >
REAL_TYPE tolerance( REAL_TYPE const a, REAL_TYPE const b, REAL_TYPE const ndigits )
{
  return std::numeric_limits< double >::epsilon() * std::max( fabs( a ), fabs( b ) ) * pow( 10, ndigits );
}

//******************************************************************************
template< int numReactions, int numSpecies >
struct ComputeReactionRatesTestData
{
  double speciesConcentration[numSpecies];
  double reactionRates[numReactions] = { 0.0 };
  CArrayWrapper< double, numReactions, numSpecies > reactionRatesDerivatives;
};

template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void computeReactionRatesTest( PARAMS_DATA const & params,
                               REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies()],
                               REAL_TYPE const (&expectedReactionRates)[PARAMS_DATA::numReactions()],
                               REAL_TYPE const (&expectedReactionRatesDerivatives)[PARAMS_DATA::numReactions()][PARAMS_DATA::numSpecies()] )
{
  using KineticReactionsType = reactionsSystems::KineticReactions< REAL_TYPE,
                                                                   int,
                                                                   int,
                                                                   LOGE_CONCENTRATION >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies();
  constexpr int numReactions = PARAMS_DATA::numReactions();

  double const temperature = 298.15;
  ComputeReactionRatesTestData< numReactions, numSpecies > data;

  REAL_TYPE magScale = 0;
  for( int r=0; r<numSpecies; ++r )
  {
    magScale = std::max( magScale, fabs( initialSpeciesConcentration[r] ));
  }


  if constexpr( LOGE_CONCENTRATION )
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      data.speciesConcentration[i] = log( initialSpeciesConcentration[i] );
    }
  }
  else
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      data.speciesConcentration[i] = initialSpeciesConcentration[i];
    }
  }


  pmpl::genericKernelWrapper( 1, &data, [params, temperature] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    KineticReactionsType::computeReactionRates( temperature,
                                                params,
                                                dataCopy->speciesConcentration,
                                                dataCopy->reactionRates,
                                                dataCopy->reactionRatesDerivatives );
  });

  for( int r=0; r<numReactions; ++r )
  {
    EXPECT_NEAR( data.reactionRates[r], expectedReactionRates[r], std::max( magScale, fabs( expectedReactionRates[r] ) ) * 1.0e-8 );
  }



  for( int r = 0; r < numReactions; ++r )
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      if constexpr( LOGE_CONCENTRATION )
      {
        data.reactionRatesDerivatives( r, i ) = data.reactionRatesDerivatives( r, i ) * exp( -data.speciesConcentration[i] );
      }
      EXPECT_NEAR( data.reactionRatesDerivatives( r, i ), expectedReactionRatesDerivatives[r][i], std::max( magScale, fabs( expectedReactionRatesDerivatives[r][i] ) ) * 1.0e-8 );
    }
  }
}


//******************************************************************************
template< int numSpecies >
struct ComputeSpeciesRatesTestData
{
  double speciesConcentration[numSpecies];
  double speciesRates[numSpecies] = { 0.0 };
  CArrayWrapper< double, numSpecies, numSpecies > speciesRatesDerivatives;
};

template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void computeSpeciesRatesTest( PARAMS_DATA const & params,
                              REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies()],
                              REAL_TYPE const (&expectedSpeciesRates)[PARAMS_DATA::numSpecies()],
                              REAL_TYPE const (&expectedSpeciesRatesDerivatives)[PARAMS_DATA::numSpecies()][PARAMS_DATA::numSpecies()] )
{

  using KineticReactionsType = reactionsSystems::KineticReactions< REAL_TYPE,
                                                                   int,
                                                                   int,
                                                                   LOGE_CONCENTRATION >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies();

  double const temperature = 298.15;
  ComputeSpeciesRatesTestData< numSpecies > data;

  if constexpr( LOGE_CONCENTRATION )
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      data.speciesConcentration[i] = log( initialSpeciesConcentration[i] );
    }
  }
  else
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      data.speciesConcentration[i] = initialSpeciesConcentration[i];
    }
  }

  pmpl::genericKernelWrapper( 1, &data, [params, temperature] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    KineticReactionsType::computeSpeciesRates( temperature,
                                               params,
                                               dataCopy->speciesConcentration,
                                               dataCopy->speciesRates,
                                               dataCopy->speciesRatesDerivatives );
  });


  for( int r=0; r<numSpecies; ++r )
  {
    EXPECT_NEAR( data.speciesRates[r], expectedSpeciesRates[r], 1.0e-8 );
  }

  for( int i = 0; i < numSpecies; ++i )
  {
    for( int j = 0; j < numSpecies; ++j )
    {
      if constexpr( LOGE_CONCENTRATION )
      {
        data.speciesRatesDerivatives( i, j ) = data.speciesRatesDerivatives( i, j ) * exp( -data.speciesConcentration[j] );
      }
      EXPECT_NEAR( data.speciesRatesDerivatives( i, j ), expectedSpeciesRatesDerivatives[i][j], 1.0e-8 );
    }
  }
}

//******************************************************************************
template< int numSpecies >
struct TimeStepTestData
{
  double speciesConcentration[numSpecies];
  double time = 0.0;
};

template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void timeStepTest( PARAMS_DATA const & params,
                   REAL_TYPE const dt,
                   int const numSteps,
                   REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies()],
                   REAL_TYPE const (&expectedSpeciesConcentrations)[PARAMS_DATA::numSpecies()] )
{
  using KineticReactionsType = reactionsSystems::KineticReactions< REAL_TYPE,
                                                                   int,
                                                                   int,
                                                                   LOGE_CONCENTRATION >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies();
  double const temperature = 298.15;
  TimeStepTestData< numSpecies > data;

  if constexpr( LOGE_CONCENTRATION )
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      data.speciesConcentration[i] = log( initialSpeciesConcentration[i] );
    }
  }
  else
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      data.speciesConcentration[i] = initialSpeciesConcentration[i];
    }
  }


  data.time = 0.0;

  pmpl::genericKernelWrapper( 1, &data, [params, temperature, dt, numSteps] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    double speciesConcentration_n[numSpecies];
    double speciesRates[numSpecies] = { 0.0 };
    CArrayWrapper< double, numSpecies, numSpecies > speciesRatesDerivatives;

    for( int t = 0; t < numSteps; ++t )
    {
      printf("Time step %d \n ", t);
      for( int i=0; i<numSpecies; ++i )
      {
        speciesConcentration_n[i] = dataCopy->speciesConcentration[i];
      }
      KineticReactionsType::timeStep( dt,
                                      temperature,
                                      params,
                                      speciesConcentration_n,
                                      dataCopy->speciesConcentration,
                                      speciesRates,
                                      speciesRatesDerivatives );
      dataCopy->time += dt;
    }
  });


  EXPECT_NEAR( data.time, dt*numSteps, 1.0e-8 );

  for( int i = 0; i < numSpecies; ++i )
  {
    if constexpr( LOGE_CONCENTRATION )
    {
      data.speciesConcentration[i] = exp( data.speciesConcentration[i] );
    }
    EXPECT_NEAR( data.speciesConcentration[i], expectedSpeciesConcentrations[i], 1.0e-4 );
  }
}

} // namespace unitTest_utilities
} // namespace hpcReact
