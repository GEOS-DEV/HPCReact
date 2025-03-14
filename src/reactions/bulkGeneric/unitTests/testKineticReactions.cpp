
#include "../KineticReactions.hpp"
#include "../ParametersPredefined.hpp"
#include "common/DirectSystemSolve.hpp"
#include "common/macros.hpp"
#include "common/printers.hpp"

#include <gtest/gtest.h>


using namespace hpcReact;
using namespace hpcReact::bulkGeneric;

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
TEST( testKineticReactions, computeReactionRatesTest_simpleTestRateParams )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedReactionRates[] = { 1.0, 0.25 };
  double const expectedReactionRatesDerivatives[][5] =
  { { 2.0, -0.5, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.5, 0.25, 0.0 } };
  computeReactionRatesTest< double, false >( simpleTestRateParams.kineticReactionsParameters(),
                                             initialSpeciesConcentration,
                                             expectedReactionRates,
                                             expectedReactionRatesDerivatives );
  computeReactionRatesTest< double, true >( simpleTestRateParams,
                                            initialSpeciesConcentration,
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
}

TEST( testKineticReactions, computeReactionRatesTest_carbonateSystem )
{
  double const initialSpeciesConcentration[18] =
  {
    1.0e-16, // OH-
    1.0e-16, // CO2
    1.0e-16, // CO3-2
    1.0e-16, // H2CO3
    1.0e-16, // CaHCO3+
    1.0e-16, // CaCO3
    1.0e-16, // CaSO4
    1.0e-16, // CaCl+
    1.0e-16, // CaCl2
    1.0e-16, // MgSO4
    1.0e-16, // NaSO4-
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89, // Cl-
    1.65e-2, // Mg+2
    1.09 // Na+1
  };

  double const expectedReactionRates[11] = { -0.001424736, -12610.7392, -0.175591624, -473.6096, -269197.19999999984, -0.00012441275624000003, -18012.914999999986, -1.56526019999999e6,
                                             -346983.07769999903, -14247.58499999999, -2.316271799999999e6 };
  double const expectedReactionRatesDerivatives[11][18] =
  {
    { 5.264e10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000014, 0, 0, 0, 0, 0, 0 },
    { 0, 0.039, 0, 0, 0, 0, 0, 0, 0, 0, 0, -33539.2, -33539.2, 0, 0, 0, 0, 0 },
    { 0, 0, 3.76e9, 0, 0, 0, 0, 0, 0, 0, 0, 1.e-6, -0.467, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0.57, 0, 0, 0, 0, 0, 0, 0, -1259.6, -1259.6, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 1.5e6, 0, 0, 0, 0, 0, 0, 0, -715950., -6.956e6, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 37600., 0, 0, 0, 0, 0, 1.e-11, -0.000330885, -0.0032148000000000003, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 100000., 0, 0, 0, 0, 0, 0, -465449.99999999994, -561150., 0, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 1.e8, 0, 0, 0, 0, 0, -4.0446e7, 0, -828180., 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 1.e7, 0, 0, 0, 0, -8.965971e6, 0, -367177.86, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 100000., 0, 0, 0, 0, -443850., 0, -863489.9999999999, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.e7, 0, 0, 0, -7.2158e7, 0, 0, -2.12502e6 }
  };

  computeReactionRatesTest< double, false >( carbonateSystem,
                                             initialSpeciesConcentration,
                                             expectedReactionRates,
                                             expectedReactionRatesDerivatives );
  computeReactionRatesTest< double, true >( carbonateSystem,
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

TEST( testKineticReactions, computeSpeciesRatesTest_simpleTestRateParams )
{
  double const initialSpeciesConcentration[5] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedSpeciesRates[5] = { -2.0, 1.0, 0.75, -0.25, 0.5 };
  double const expectedSpeciesRatesDerivatives[5][5] = { { -4.0, 1.0, 0.0, 0.0, 0.0 },
    {  2.0, -0.5, 0.0, 0.0, 0.0 },
    {  2.0, -0.5, -0.5, -0.25, 0.0 },
    {  0.0, 0.0, -0.5, -0.25, 0.0 },
    {  0.0, 0.0, 1.0, 0.5, 0.0 } };

  computeSpeciesRatesTest< double, false >( simpleTestRateParams,
                                            initialSpeciesConcentration,
                                            expectedSpeciesRates,
                                            expectedSpeciesRatesDerivatives );

  computeSpeciesRatesTest< double, true >( simpleTestRateParams,
                                           initialSpeciesConcentration,
                                           expectedSpeciesRates,
                                           expectedSpeciesRatesDerivatives );

}

// TEST( testKineticReactions, computeSpeciesRatesTest_carbonateSystem )
// {
//   double const initialSpeciesConcentration[18] =
//   {
//     1.0e-16, // OH-
//     1.0e-16, // CO2
//     1.0e-16, // CO3-2
//     1.0e-16, // H2CO3
//     1.0e-16, // CaHCO3+
//     1.0e-16, // CaCO3
//     1.0e-16, // CaSO4
//     1.0e-16, // CaCl+
//     1.0e-16, // CaCl2
//     1.0e-16, // MgSO4
//     1.0e-16, // NaSO4-
//     3.76e-1, // H+
//     3.76e-1, // HCO3-
//     3.87e-2, // Ca+2
//     3.21e-2, // SO4-2
//     1.89, // Cl-
//     1.65e-2, // Mg+2
//     1.09 // Na+1
//   };

//   double const expectedSpeciesRates[18] = { 0 };
//   double const expectedSpeciesRatesDerivatives[18][18] = {{ 0}};

//   computeSpeciesRatesTest< double, false >( carbonateSystem,
//                                             initialSpeciesConcentration,
//                                             expectedSpeciesRates,
//                                             expectedSpeciesRatesDerivatives );


// }

//******************************************************************************
template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void timeStepTest( PARAMS_DATA const & params,
                   REAL_TYPE const dt,
                   int const numSteps,
                   REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies],
                   REAL_TYPE const (&expectedSpeciesConcentrations)[PARAMS_DATA::numSpecies] )
{
  using KineticReactionsType = KineticReactions< double,
                                                 int,
                                                 int,
                                                 LOGE_CONCENTRATION >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies;
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

TEST( testKineticReactions, testTimeStep )
{
  double const initialSpeciesConcentration[5] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedSpeciesConcentrations[5] = { 3.92138293924124e-01, 3.03930853037938e-01, 5.05945480771998e-01, 7.02014627734060e-01, 5.95970744531880e-01 };

  timeStepTest< double, false >( simpleTestRateParams,
                                 2.0,
                                 10,
                                 initialSpeciesConcentration,
                                 expectedSpeciesConcentrations );

  // ln(c) as the primary variable results in a singular system.
  // timeStepTest< double, true >( simpleTestRateParams,
  //                               2.0,
  //                               10,
  //                               initialSpeciesConcentration,
  //                               expectedSpeciesConcentrations );
}


TEST( testKineticReactions, testTimeStep_carbonateSystem )
{
  double const initialSpeciesConcentration[18] =
  {
    1.0e-16, // OH-
    1.0e-16, // CO2
    1.0e-16, // CO3-2
    1.0e-16, // H2CO3
    1.0e-16, // CaHCO3+
    1.0e-16, // CaCO3
    1.0e-16, // CaSO4
    1.0e-16, // CaCl+
    1.0e-16, // CaCl2
    1.0e-16, // MgSO4
    1.0e-16, // NaSO4-
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89, // Cl-
    1.65e-2, // Mg+2
    1.09 // Na+1
  };

  double const expectedSpeciesConcentrations[18] =
  { 2.327841695586879e-11, // OH-
    3.745973700632716e-01, // CO2
    3.956656978189456e-11, // CO3-2
    9.629355924567627e-04, // H2CO3
    6.739226982791492e-05, // CaHCO3+
    1.065032288527957e-09, // CaCO3
    5.298329882666738e-03, // CaSO4
    5.844517547638333e-03, // CaCl+
    1.277319392670652e-02, // CaCl2
    6.618125707964991e-03, // MgSO4
    1.769217213462983e-02, // NaSO4-
    4.396954721488358e-04, // H+
    3.723009698453808e-04, // HCO3-
    1.471656530812871e-02, // Ca+2
    2.491372274738741e-03, // SO4-2
    1.858609094598949e+00, // Cl-
    9.881874292035110e-03, // Mg+2
    1.072307827865370e+00 // Na+1
  };

  timeStepTest< double, false >( carbonateSystem,
                                 2.0,
                                 100000,
                                 initialSpeciesConcentration,
                                 expectedSpeciesConcentrations );

  // ln(c) as the primary variable results in a singular system.
  // timeStepTest< double, true >( simpleTestRateParams,
  //                               2.0,
  //                               10,
  //                               initialSpeciesConcentration,
  //                               expectedSpeciesConcentrations );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
