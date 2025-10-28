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


TEST( testKineticReactions, computeReactionRatesTest_carbonateSystemAllKinetic )
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

  double const expectedReactionRates[10] = { -0.001424736, //             OH- + H+ = H2O
                                             -12610.7392, //              CO2 + H2O = H+ + HCO3-
                                             -0.175591624, //             CO3-2 + H+ = HCO3-
                                             -269197.19999999984, //      CaHCO3+ = Ca+2 + HCO3-
                                             -18012.914999999986, //      CaSO4 = Ca+2 + SO4-2
                                             -1.56526019999999e6, //      CaCl+ = Ca+2 + Cl-
                                             -346983.07769999903, //      CaCl2 = Ca+2 + 2Cl-
                                             -14247.58499999999, //       MgSO4 = Mg+2 + SO4-2
                                             -2.316271799999999e6, //     NaSO4- = Na+ + SO4-2
                                             -4.3653599999994173e-10 // CaCO3 + H+ = Ca+2 + HCO3- (kinetic)
  };
  double const expectedReactionRatesDerivatives[10][17] =
  {
    { 5.264e10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.000014, 0, 0, 0, 0, 0, 0 },
    { 0, 0.039, 0, 0, 0, 0, 0, 0, 0, 0, -33539.2, -33539.2, 0, 0, 0, 0, 0 },
    { 0, 0, 3.76e9, 0, 0, 0, 0, 0, 0, 0, 1.e-6, -0.467, 0, 0, 0, 0, 0 },
    { 0, 0, 0, 1.5e6, 0, 0, 0, 0, 0, 0, 0, -715950., -6.956e6, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 100000., 0, 0, 0, 0, 0, 0, 0, -465449.99999999994, -561150., 0, 0, 0 },
    { 0, 0, 0, 0, 0, 1.e8, 0, 0, 0, 0, 0, 0, -4.0446e7, 0, -828180., 0, 0 },
    { 0, 0, 0, 0, 0, 0, 1.e7, 0, 0, 0, 0, 0, -8.965971e6, 0, -367177.86, 0, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 100000., 0, 0, 0, 0, 0, -443850., 0, -863489.9999999999, 0 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 1.e7, 0, 0, 0, 0, -7.2158e7, 0, 0, -2.12502e6 },
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 5.8279999999999995e-07, 1.e-11, -1.1609999999999998e-09, -1.1280000000000004e-08, 0, 0, 0, 0 }

  };

  computeReactionRatesTest< double, false >( carbonateSystemAllKinetic.kineticReactionsParameters(),
                                             initialSpeciesConcentration,
                                             expectedReactionRates,
                                             expectedReactionRatesDerivatives );
  computeReactionRatesTest< double, true >( carbonateSystemAllKinetic.kineticReactionsParameters(),
                                            initialSpeciesConcentration,
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
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


TEST( testKineticReactions, testTimeStep_carbonateSystemAllKinetic )
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

  double const expectedSpeciesConcentrations[17] =
  { 2.327841695586879e-11, // OH-
    0.37555955033916549, // CO2
    3.956656978189456e-11, // CO3-2
    6.739226982791492e-05, // CaHCO3+
    5.298329882666738e-03, // CaSO4
    5.844517547638333e-03, // CaCl+
    1.277319392670652e-02, // CaCl2
    6.618125707964991e-03, // MgSO4
    1.769217213462983e-02, // NaSO4-
    1.065032288527957e-09, // CaCO3
    4.396954721488358e-04, // H+
    3.723009698453808e-04, // HCO3-
    1.471656530812871e-02, // Ca+2
    2.491372274738741e-03, // SO4-2
    1.858609094598949e+00, // Cl-
    9.881874292035110e-03, // Mg+2
    1.072307827865370e+00 // Na+1
  };

  timeStepTest< double, false >( carbonateSystemAllKinetic.kineticReactionsParameters(),
                                 10.0,
                                 10000,
                                 initialSpeciesConcentration,
                                 expectedSpeciesConcentrations );

  // ln(c) as the primary variable results in a singular system.
  // timeStepTest< double, true >( simpleKineticTestRateParams,
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
