
#include "../KineticReactions_impl.hpp"
#include "../ParametersPredefined.hpp"
#include "common/DirectSystemSolve.hpp"
#include "common/macros.hpp"

#include <gtest/gtest.h>


using namespace hpcReact;
using namespace hpcReact::bulkGeneric;


//******************************************************************************
template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void computeReactionRatesTest( PARAMS_DATA const & params,
                               REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies],
                               REAL_TYPE const (&expectedReactionRates)[PARAMS_DATA::numReactions],
                               REAL_TYPE const (&expectedReactionRatesDerivatives)[PARAMS_DATA::numReactions][PARAMS_DATA::numSpecies] )
{
  using KineticReactionsType = KineticReactions< REAL_TYPE, 
                                                 int, 
                                                 int,
                                                 LOGE_CONCENTRATION >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies;
  constexpr int numReactions = PARAMS_DATA::numReactions;

  double const temperature = 298.15;
  double speciesConcentration[numSpecies];

  if constexpr ( LOGE_CONCENTRATION )
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
    EXPECT_NEAR( reactionRates[r], expectedReactionRates[r], 1.0e-8 );
  }

  for( int r = 0; r < numReactions; ++r )
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      if constexpr ( LOGE_CONCENTRATION )
      {
        reactionRatesDerivatives( r, i ) = reactionRatesDerivatives( r, i ) * exp( -speciesConcentration[i] );
      }
      EXPECT_NEAR( reactionRatesDerivatives( r, i ), expectedReactionRatesDerivatives[r][i], 1.0e-8 );
    }
  }
}


//******************************************************************************
TEST( bulkGeneric, computeReactionRatesTest )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedReactionRates[] = { 1.0, 0.25 };
  double const expectedReactionRatesDerivatives[][5] = { { 2.0, -0.5, 0.0, 0.0, 0.0 },
                                                          { 0.0, 0.0, 0.5, 0.25, 0.0 } };
  computeReactionRatesTest< double, false >( simpleTestRateParams,
                                             initialSpeciesConcentration,
                                             expectedReactionRates,
                                             expectedReactionRatesDerivatives );
  computeReactionRatesTest< double, true >( simpleTestRateParams,
                                            initialSpeciesConcentration,
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
}





//******************************************************************************
template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void computeSpeciesRatesTest( PARAMS_DATA const & params,
                               REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies],
                               REAL_TYPE const (&expectedSpeciesRates)[PARAMS_DATA::numSpecies],
                               REAL_TYPE const (&expectedSpeciesRatesDerivatives)[PARAMS_DATA::numSpecies][PARAMS_DATA::numSpecies] )
{

  using KineticReactionsType = KineticReactions< REAL_TYPE, 
                                                 int, 
                                                 int,
                                                 LOGE_CONCENTRATION >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies;

  double const temperature = 298.15;
  double speciesConcentration[numSpecies];
  double speciesRates[numSpecies] = { 0.0 };
  CArrayWrapper< double, numSpecies, numSpecies > speciesRatesDerivatives;
  
  if constexpr ( LOGE_CONCENTRATION )
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
      if constexpr ( LOGE_CONCENTRATION )
      {
        speciesRatesDerivatives( i, j ) = speciesRatesDerivatives( i, j ) * exp( -speciesConcentration[j] );
      }
      EXPECT_NEAR( speciesRatesDerivatives( i, j ), expectedSpeciesRatesDerivatives[i][j], 1.0e-8 );
    }
  }
}

TEST( bulkGeneric, computeSpeciesRatesTest )
{
  double const initialSpeciesConcentration[5] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedSpeciesRates[5] = { -2.0, 1.0, 0.75, -0.25, 0.5 };
  double const expectedSpeciesRatesDerivatives[5][5] = { { -4.0,  1.0,  0.0,   0.0, 0.0 }, 
                                                         {  2.0, -0.5,  0.0,   0.0, 0.0 }, 
                                                         {  2.0, -0.5, -0.5, -0.25, 0.0 }, 
                                                         {  0.0,  0.0, -0.5, -0.25, 0.0 }, 
                                                         {  0.0,  0.0,  1.0,   0.5, 0.0 } };

  computeSpeciesRatesTest< double, false >( simpleTestRateParams,
                                             initialSpeciesConcentration,
                                             expectedSpeciesRates,
                                             expectedSpeciesRatesDerivatives );

  computeSpeciesRatesTest< double, true >( simpleTestRateParams,
                                              initialSpeciesConcentration,
                                              expectedSpeciesRates,
                                              expectedSpeciesRatesDerivatives );
 
}

//******************************************************************************
template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void timeStepTest( PARAMS_DATA const & params,
                   REAL_TYPE const dt,
                   int const numSteps,
                   REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies],
                   REAL_TYPE const (&expectedSpeciesConcentrations)[PARAMS_DATA::numSpecies],
                   REAL_TYPE const (&expectedSpeciesRates)[PARAMS_DATA::numSpecies],
                   REAL_TYPE const (&expectedSpeciesRatesDerivatives)[PARAMS_DATA::numSpecies][PARAMS_DATA::numSpecies] )
{
  using KineticReactionsType = KineticReactions< double, 
                                                 int, 
                                                 int,
                                                 LOGE_CONCENTRATION >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies;
  double const temperature = 298.15;

  double speciesConcentration[numSpecies];
  if constexpr ( LOGE_CONCENTRATION )
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

  double speciesRates[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

  double time = 0.0;
  for( int t = 0; t < numSteps; ++t )
  {
    CArrayWrapper<double,numSpecies,numSpecies> speciesRateDerivatives;
    double speciesConcentration_n[] = { speciesConcentration[0], 
                                               speciesConcentration[1], 
                                               speciesConcentration[2], 
                                               speciesConcentration[3], 
                                               speciesConcentration[4] };



    KineticReactionsType::timeStep( dt, 
                                    temperature, 
                                    params, 
                                    speciesConcentration_n, 
                                    speciesConcentration, 
                                    speciesRates,
                                    speciesRateDerivatives );

    time += dt;
  }

  HPCREACT_UNUSED_VAR( expectedSpeciesRates );
  HPCREACT_UNUSED_VAR( expectedSpeciesConcentrations);
  HPCREACT_UNUSED_VAR( expectedSpeciesRatesDerivatives );
  EXPECT_NEAR( time, dt*numSteps, 1.0e-8 );

  for( int i = 0; i < numSpecies; ++i )
  {
    if constexpr ( LOGE_CONCENTRATION )
    {
      speciesConcentration[i] = exp( speciesConcentration[i] );
    }
    EXPECT_NEAR( speciesConcentration[i], expectedSpeciesConcentrations[i], 1.0e-4 );

    // if constexpr ( LOGE_CONCENTRATION )
    // {
    //   speciesRatesDerivatives( i, j ) = speciesRatesDerivatives( i, j ) * exp( -speciesConcentration[j] );
    // }
    // EXPECT_NEAR( speciesRatesDerivatives( i, j ), expectedSpeciesRatesDerivatives[i][j], 1.0e-8 );

    
  }
}

TEST( bulkGeneric, testTimeStep )
{
  double const initialSpeciesConcentration[5] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedSpeciesConcentrations[5] = { 3.92138294e-01, 3.03930853e-01, 5.05945481e-01, 7.02014628e-01, 5.95970745e-01 };
  double const expectedSpeciesRates[5] = { -2.0, 1.0, 0.75, -0.25, 0.5 };
  double const expectedSpeciesRatesDerivatives[5][5] = { { -4.0,  1.0,  0.0,   0.0, 0.0 }, 
                                                         {  2.0, -0.5,  0.0,   0.0, 0.0 }, 
                                                         {  2.0, -0.5, -0.5, -0.25, 0.0 }, 
                                                         {  0.0,  0.0, -0.5, -0.25, 0.0 }, 
                                                         {  0.0,  0.0,  1.0,   0.5, 0.0 } };

  timeStepTest< double, false >( simpleTestRateParams,
                                 2.0,
                                 10,
                                 initialSpeciesConcentration,
                                 expectedSpeciesConcentrations,
                                 expectedSpeciesRates,
                                 expectedSpeciesRatesDerivatives );

  // ln(c) as the primary variable results in a singular system.
  // timeStepTest< double, true >( simpleTestRateParams,
  //                               2.0,
  //                               10,
  //                               initialSpeciesConcentration,
  //                               expectedSpeciesConcentrations,
  //                               expectedSpeciesRates,
  //                               expectedSpeciesRatesDerivatives );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
