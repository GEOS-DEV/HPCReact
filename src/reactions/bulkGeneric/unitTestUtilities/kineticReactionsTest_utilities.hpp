#pragma once
#include "common/DirectSystemSolve.hpp"
#include "common/macros.hpp"
#include "common/printers.hpp"
#include "common/CArrayWrapper.hpp"
#include "../KineticReactions.hpp"

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
template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void computeReactionRatesTest( PARAMS_DATA const & params,
                               REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies()],
                               REAL_TYPE const (&expectedReactionRates)[PARAMS_DATA::numReactions()],
                               REAL_TYPE const (&expectedReactionRatesDerivatives)[PARAMS_DATA::numReactions()][PARAMS_DATA::numSpecies()] )
{
  using KineticReactionsType = bulkGeneric::KineticReactions< REAL_TYPE,
                                                              int,
                                                              int,
                                                              LOGE_CONCENTRATION >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies();
  constexpr int numReactions = PARAMS_DATA::numReactions();

  double const temperature = 298.15;
  double speciesConcentration[numSpecies];

  REAL_TYPE magScale = 0;
  for( int r=0; r<numSpecies; ++r )
  {
    magScale = std::max( magScale, fabs( initialSpeciesConcentration[r] ));
  }


  if constexpr( LOGE_CONCENTRATION )
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      speciesConcentration[i] = log( initialSpeciesConcentration[i] );
    }
  }
  else
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      speciesConcentration[i] = initialSpeciesConcentration[i];
    }
  }

  double reactionRates[numReactions] = { 0.0 };
  CArrayWrapper< double, numReactions, numSpecies > reactionRatesDerivatives;


  KineticReactionsType::computeReactionRates( temperature,
                                              params,
                                              speciesConcentration,
                                              reactionRates,
                                              reactionRatesDerivatives );



  for( int r=0; r<numReactions; ++r )
  {
    EXPECT_NEAR( reactionRates[r], expectedReactionRates[r], std::max( magScale, fabs( expectedReactionRates[r] ) ) * 1.0e-8 );
  }



  for( int r = 0; r < numReactions; ++r )
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      if constexpr( LOGE_CONCENTRATION )
      {
        reactionRatesDerivatives( r, i ) = reactionRatesDerivatives( r, i ) * exp( -speciesConcentration[i] );
      }
      EXPECT_NEAR( reactionRatesDerivatives( r, i ), expectedReactionRatesDerivatives[r][i], std::max( magScale, fabs( expectedReactionRatesDerivatives[r][i] ) ) * 1.0e-8 );
    }
  }
}


//******************************************************************************
template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void computeSpeciesRatesTest( PARAMS_DATA const & params,
                              REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies()],
                              REAL_TYPE const (&expectedSpeciesRates)[PARAMS_DATA::numSpecies()],
                              REAL_TYPE const (&expectedSpeciesRatesDerivatives)[PARAMS_DATA::numSpecies()][PARAMS_DATA::numSpecies()] )
{

  using KineticReactionsType = bulkGeneric::KineticReactions< REAL_TYPE,
                                                              int,
                                                              int,
                                                              LOGE_CONCENTRATION >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies();

  double const temperature = 298.15;
  double speciesConcentration[numSpecies];
  double speciesRates[numSpecies] = { 0.0 };
  CArrayWrapper< double, numSpecies, numSpecies > speciesRatesDerivatives;

  if constexpr( LOGE_CONCENTRATION )
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      speciesConcentration[i] = log( initialSpeciesConcentration[i] );
    }
  }
  else
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      speciesConcentration[i] = initialSpeciesConcentration[i];
    }
  }

  KineticReactionsType::computeSpeciesRates( temperature,
                                             params,
                                             speciesConcentration,
                                             speciesRates,
                                             speciesRatesDerivatives );


  for( int r=0; r<numSpecies; ++r )
  {
    EXPECT_NEAR( speciesRates[r], expectedSpeciesRates[r], 1.0e-8 );
  }

  for( int i = 0; i < numSpecies; ++i )
  {
    for( int j = 0; j < numSpecies; ++j )
    {
      if constexpr( LOGE_CONCENTRATION )
      {
        speciesRatesDerivatives( i, j ) = speciesRatesDerivatives( i, j ) * exp( -speciesConcentration[j] );
      }
      EXPECT_NEAR( speciesRatesDerivatives( i, j ), expectedSpeciesRatesDerivatives[i][j], 1.0e-8 );
    }
  }
}

//******************************************************************************
template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void timeStepTest( PARAMS_DATA const & params,
                   REAL_TYPE const dt,
                   int const numSteps,
                   REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies()],
                   REAL_TYPE const (&expectedSpeciesConcentrations)[PARAMS_DATA::numSpecies()] )
{
  using KineticReactionsType = bulkGeneric::KineticReactions< REAL_TYPE,
                                                              int,
                                                              int,
                                                              LOGE_CONCENTRATION >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies();
  double const temperature = 298.15;

  double speciesConcentration[numSpecies];
  if constexpr( LOGE_CONCENTRATION )
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      speciesConcentration[i] = log( initialSpeciesConcentration[i] );
    }
  }
  else
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      speciesConcentration[i] = initialSpeciesConcentration[i];
    }
  }


  double time = 0.0;
  for( int t = 0; t < numSteps; ++t )
  {
    double speciesConcentration_n[numSpecies];
    for( int i=0; i<numSpecies; ++i )
    {
      speciesConcentration_n[i] = speciesConcentration[i];
    }

    double speciesRates[numSpecies] = { 0.0 };
    CArrayWrapper< double, numSpecies, numSpecies > speciesRatesDerivatives;

    KineticReactionsType::timeStep( dt,
                                    temperature,
                                    params,
                                    speciesConcentration_n,
                                    speciesConcentration,
                                    speciesRates,
                                    speciesRatesDerivatives );

    time += dt;
  }

  EXPECT_NEAR( time, dt*numSteps, 1.0e-8 );

  for( int i = 0; i < numSpecies; ++i )
  {
    if constexpr( LOGE_CONCENTRATION )
    {
      speciesConcentration[i] = exp( speciesConcentration[i] );
    }
    EXPECT_NEAR( speciesConcentration[i], expectedSpeciesConcentrations[i], 1.0e-4 );
  }
}

} // namespace unitTest_utilities
} // namespace hpcReact