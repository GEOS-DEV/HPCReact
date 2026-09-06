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


#include "reactions/unitTestUtilities/kineticReactionsTestUtilities.hpp"
#include "../GeochemicalSystems.hpp"
#include <gtest/gtest.h>


using namespace hpcReact;
using namespace hpcReact::geochemistry;
using namespace hpcReact::unitTest_utilities;


TEST( testKineticReactions, computeReactionRatesTest_carbonateSystemAllKinetic_Identity )
{
  double const initialSpeciesConcentration[17] =
  {
    1.0e-16, // OH-
    1.0e-16, // CO2
    1.0e-16, // CO3-2
    1.0e-16, // CaHCO3+
    1.0e-16, // CaSO4
    1.0e-16, // CaCl+
    1.0e-16, // CaCl2
    1.0e-16, // MgSO4
    1.0e-16, // NaSO4-
    1.0e-16, // CaCO3
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89, // Cl-
    1.65e-2, // Mg+2
    1.09 // Na+1
  };

  double const surfaceArea[10] = { 0.0, // OH- + H+ = H2O
                                   0.0, // CO2 + H2O = H+ + HCO3-
                                   0.0, // CO3-2 + H+ = HCO3-
                                   0.0, // CaHCO3+ = Ca+2 + HCO3-
                                   0.0, // CaSO4 = Ca+2 + SO4-2
                                   0.0, // CaCl+ = Ca+2 + Cl-
                                   0.0, // CaCl2 = Ca+2 + 2Cl-
                                   0.0, // MgSO4 = Mg+2 + SO4-2
                                   0.0, // NaSO4- = Na+ + SO4-2
                                   0.0, // CaCO3 + H+ = Ca+2 + HCO3- (kinetic)
  };

  double const expectedReactionRates[10] = { -0.001410616, //             OH- + H+ = H2O
                                             -12193.835513600001, //              CO2 + H2O = H+ + HCO3-
                                             -0.17635490400000001, //             CO3-2 + H+ = HCO3-
                                             -243047.23847999985, //      CaHCO3+ = Ca+2 + HCO3-
                                             -16044.165503999988, //      CaSO4 = Ca+2 + SO4-2
                                             -1474255.67939999, //      CaCl+ = Ca+2 + Cl-
                                             -314076.36382919899, //      CaCl2 = Ca+2 + 2Cl-
                                             -13667.512319999991, //       MgSO4 = Mg+2 + SO4-2
                                             -2311702.236599999, //     NaSO4- = Na+ + SO4-2
                                             -3.1954435199994173e-06 // CaCO3 + H+ = Ca+2 + HCO3- (kinetic)
  };
  double const expectedReactionRatesDerivatives[10][17] =
  {
    { 52640000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.4e-05, 0, 0, 0, 0, 0, 0 },
    { 0, 0.039, 0, 0, 0, 0, 0, 0, 0, 0, -32430.413600000003, -32430.413600000003, 0, 0, 0, 0, 0 },
    { 0, 0, 3760000000, 0, 0, 0, 0, 0, 0, 0, 9.9999999999999995e-07, -0.46903, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 1500000, 0, 0, 0, 0, 0, 0, 0, -646402.22999999998, -6280290.4000000004, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 100000, 0, 0, 0, 0, 0, 0, 0, -414577.91999999998, -499818.23999999999, 0, 0, 0 },
    { 0, 0, 0, 0, 0, 100000000, 0, 0, 0, 0, 0, 0, -38094462, 0, -780029.45999999996, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 10000000, 0, 0, 0, 0, 0, -8115668.3159999996, 0, -332355.94056000002, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 100000, 0, 0, 0, 0, 0, -425779.20000000001, 0, -828334.07999999996, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 10000000, 0, 0, 0, 0, -72015646, 0, 0, -2120827.7399999998 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.8279999999999995e-03, 1.5500000000000001e-18, -8.4985199999999993e-06, -8.2569599999999996e-05, 0, 0, 0, 0 }

  };

  using ActivityType = carbonateIdentityActivityType;


  computeReactionRatesTest< double,
                            false,
                            ActivityType >( carbonateSystemAllKinetic.kineticReactionsParameters(),
                                            carbonateIdentityActivityParams,
                                            initialSpeciesConcentration,
                                            surfaceArea, // No use. Just to pass something here
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
  computeReactionRatesTest< double,
                            true,
                            ActivityType >( carbonateSystemAllKinetic.kineticReactionsParameters(),
                                            carbonateIdentityActivityParams,
                                            initialSpeciesConcentration,
                                            surfaceArea, // No use. Just to pass something here
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
}

TEST( testKineticReactions, computeReactionRatesQuotientTest_carbonateSystem_Identity )
{
  double const initialSpeciesConcentration[16] =
  {
    1.0e-16, // OH-
    1.0e-16, // CO2
    1.0e-16, // CO3-2
    1.0e-16, // CaHCO3+
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

  double const surfaceArea[1] = { 1e2 }; // CaCO3 + H+ = Ca+2 + HCO3- (kinetic)

  double const expectedReactionRates[1] = { 1.5491501480000001 }; // CaCO3 + H+ = Ca+2 + HCO3- (kinetic)

  double const expectedReactionRatesDerivatives[1][16] =
  {
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0022602446808510633, -0.0022602446808510633, -0.02196, 0, 0, 0, 0 }
  };

  using ActivityType = carbonateNosolidIdentityActivityType;

  computeReactionRatesTest< double,
                            false,
                            ActivityType >( carbonateSystem.kineticReactionsParameters(),
                                            carbonateNosolidIdentityActivityParams,
                                            initialSpeciesConcentration,
                                            surfaceArea,
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
  computeReactionRatesTest< double,
                            true,
                            ActivityType >( carbonateSystem.kineticReactionsParameters(),
                                            carbonateNosolidIdentityActivityParams,
                                            initialSpeciesConcentration,
                                            surfaceArea,
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
}

//******************************************************************************

/**
 * @brief Compute the calcite reaction rate for a given set of concentrations.
 */
template< bool LOGE_CONCENTRATION >
double calciteReactionRate( double const (&speciesConcentration)[16],
                            double const surfaceAreaValue )
{
  using ActivityType = carbonateNosolidActivityType;
  using KineticReactionsType = reactionsSystems::KineticReactions< double,
                                                                   int,
                                                                   int,
                                                                   ActivityType,
                                                                   LOGE_CONCENTRATION >;

  auto const params = carbonateSystem.kineticReactionsParameters();

  // Captured by value: a namespace-scope constexpr is host-only inside a device lambda.
  auto const activityParams = hpcReact::geochemistry::carbonateNosolidActivityParamsEQ36;

  ComputeReactionRatesTestData< 1, 16 > data;
  for( int i = 0; i < 16; ++i )
  {
    data.speciesConcentration[i] = LOGE_CONCENTRATION ? log( speciesConcentration[i] )
                                                      : speciesConcentration[i];
  }
  data.surfaceArea[0] = surfaceAreaValue;

  pmpl::genericKernelWrapper( 1, &data, [params, activityParams] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    KineticReactionsType::computeReactionRates( 298.15,
                                                params,
                                                activityParams,
                                                dataCopy->speciesConcentration,
                                                dataCopy->surfaceArea,
                                                dataCopy->reactionRates,
                                                dataCopy->reactionRatesDerivatives );
  } );

  return data.reactionRates[0];
}

// The rate of the calcite reaction, verified against EQ3NR.
//
// The brine is the converged state of eq36Database/carbonate.3o, for which EQ3NR reports a calcite
// saturation state of log Q/K = -4.14608. That fixes the expected rate through r = k * A * (1 - Q/K)
// and so compares this code's activity model, ion activity product and equilibrium constant against
// EQ3NR's own saturation calculation.
TEST( testKineticReactions, computeReactionRatesVsEQ36_carbonateSystem_Bdot )
{
  // EQ3NR converged molalities, in this code's species order.
  double const speciesConcentration[16] =
  {
    2.6702e-11, // OH-
    3.7534e-01, // CO2
    2.3166e-10, // CO3-2
    4.7225e-05, // CaHCO3+
    1.6821e-03, // CaSO4
    2.3750e-03, // CaCl+
    2.0222e-03, // CaCl2
    2.0650e-03, // MgSO4
    1.3357e-02, // NaSO4-
    6.5867e-04, // H+
    6.1144e-04, // HCO3-
    3.2573e-02, // Ca+2
    1.4996e-02, // SO4-2
    1.8836e+00, // Cl-
    1.4435e-02, // Mg+2
    1.0766e+00 // Na+1
  };

  double const surfaceArea = 1.0;

  // EQ3NR 'Calcite  -4.14608', from the saturation states of the pure solids.
  double const eq36Log10QOverK = -4.14608;

  double const expectedReactionRate =
    carbonateSystem.kineticReactionsParameters().rateConstantForward( 0 ) * surfaceArea *
    ( 1.0 - pow( 10.0, eq36Log10QOverK ) );

  EXPECT_NEAR( calciteReactionRate< false >( speciesConcentration, surfaceArea ),
               expectedReactionRate,
               1.0e-7 * expectedReactionRate );

  EXPECT_NEAR( calciteReactionRate< true >( speciesConcentration, surfaceArea ),
               expectedReactionRate,
               1.0e-7 * expectedReactionRate );
}

//******************************************************************************


// TEST( testKineticReactions, computeSpeciesRatesTest_carbonateSystemAllKinetic )
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

//   computeSpeciesRatesTest< double, false >( carbonateSystemAllKinetic,
//                                             initialSpeciesConcentration,
//                                             expectedSpeciesRates,
//                                             expectedSpeciesRatesDerivatives );


// }


// TEST( testKineticReactions, testTimeStep_carbonateSystemAllKinetic )
// {
//   double const initialSpeciesConcentration[17] =
//   {
//     1.0e-16, // OH-
//     1.0e-16, // CO2
//     1.0e-16, // CO3-2
//     1.0e-16, // CaHCO3+
//     1.0e-16, // CaSO4
//     1.0e-16, // CaCl+
//     1.0e-16, // CaCl2
//     1.0e-16, // MgSO4
//     1.0e-16, // NaSO4-
//     1.0e-16, // CaCO3
//     3.76e-1, // H+
//     3.76e-1, // HCO3-
//     3.87e-2, // Ca+2
//     3.21e-2, // SO4-2
//     1.89, // Cl-
//     1.65e-2, // Mg+2
//     1.09 // Na+1
//   };

//   double const expectedSpeciesConcentrations[17] =
//   { 2.327841695586879e-11, // OH-
//     0.37555955033916549, // CO2
//     3.956656978189456e-11, // CO3-2
//     6.739226982791492e-05, // CaHCO3+
//     5.298329882666738e-03, // CaSO4
//     5.844517547638333e-03, // CaCl+
//     1.277319392670652e-02, // CaCl2
//     6.618125707964991e-03, // MgSO4
//     1.769217213462983e-02, // NaSO4-
//     1.065032288527957e-09, // CaCO3
//     4.396954721488358e-04, // H+
//     3.723009698453808e-04, // HCO3-
//     1.471656530812871e-02, // Ca+2
//     2.491372274738741e-03, // SO4-2
//     1.858609094598949e+00, // Cl-
//     9.881874292035110e-03, // Mg+2
//     1.072307827865370e+00 // Na+1
//   };

//   using ActivityType = carbonateIdentityActivityType;

//   timeStepTest< double,
//                 false,
//                 ActivityType >( carbonateSystemAllKinetic.kineticReactionsParameters(),
//                                 ActivityType::Params(),
//                                  10.0,
//                                  10000,
//                                  initialSpeciesConcentration,
//                                  expectedSpeciesConcentrations );

// ln(c) as the primary variable results in a singular system.
// timeStepTest< double, true >( simpleKineticTestRateParams,
//                               2.0,
//                               10,
//                               initialSpeciesConcentration,
//                               expectedSpeciesConcentrations );
//}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
