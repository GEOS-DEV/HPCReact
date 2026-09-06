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


#include "reactions/unitTestUtilities/equilibriumReactionsTestUtilities.hpp"
#include "../GeochemicalSystems.hpp"

using namespace hpcReact;
using namespace hpcReact::geochemistry;
using namespace hpcReact::unitTest_utilities;


// //******************************************************************************
// TEST( testEquilibriumReactions, testEnforceEquilibrium )
// {
//   double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
//   double const expectedSpeciesConcentrations[5] = { 3.92138294e-01, 3.03930853e-01, 5.05945481e-01, 7.02014628e-01, 5.95970745e-01 };


//   std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
//   testEnforceEquilibrium< double, 2 >( simpleTestRateParams.equilibriumReactionsParameters(),
//                                        initialSpeciesConcentration,
//                                        expectedSpeciesConcentrations );

// }


//******************************************************************************
TEST( testEquilibriumReactions, testcarbonateSystemAllEquilibrium_Identity )
{
  using namespace hpcReact::geochemistry;

  double const initialSpeciesConcentration[17] =
  {
    1.0e-16, // OH-
    1.0e-16, // CO2
    1.0e-16, // CO3-2
    //1.0e-16, // H2CO3
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
  { 1.7631991300262666e-11, // OH-
    0.37553965856049809, // CO2
    2.4208013881700254e-11, // CO3-2
    5.1052345859666575e-05, // CaHCO3+
    0.0050728241823463005, // CaSO4
    0.0058054702754571424, // CaCl+
    0.012170717821098593, // CaCl2
    0.0065270307598572714, // MgSO4
    0.017963901804350708, // NaSO4-
    0.00011324449476863112, // CaCO3
    0.00057358597611036207, // H+
    0.00029604457466600952, // HCO3-
    0.015486690880470161, // Ca+2
    0.0025362432534460182, // SO4-2
    1.8598530940823459, // Cl-
    0.0099729692401428292, // Mg+2
    1.0720360981956494 // Na+1
  };

  std::cout<<" RESIDUAL_FORM 0:"<<std::endl;
  testEnforceEquilibrium< double, 0, carbonateIdentityActivityType >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                                      hpcReact::geochemistry::carbonateIdentityActivityParams,
                                                                      initialSpeciesConcentration,
                                                                      expectedSpeciesConcentrations );

  // std::cout<<" RESIDUAL_FORM 1:"<<std::endl;
  // testEnforceEquilibrium< double, 1 >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
  //                                      initialSpeciesConcentration,
  //                                      expectedSpeciesConcentrations );

  std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
  testEnforceEquilibrium< double, 2, carbonateIdentityActivityType >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                                      hpcReact::geochemistry::carbonateIdentityActivityParams,
                                                                      initialSpeciesConcentration,
                                                                      expectedSpeciesConcentrations );

}


TEST( testEquilibriumReactions, testcarbonateSystemAllEquilibrium2_Identity )
{


  static constexpr int numPrimarySpecies = hpcReact::geochemistry::carbonateSystemAllEquilibrium.numPrimarySpecies();
  static constexpr int numSpecies = hpcReact::geochemistry::carbonateSystemAllEquilibrium.numSpecies();

  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< double,
                                                                           int,
                                                                           int,
                                                                           Identity< double, int, SpeciatedIonicStrength< double, int, numSpecies > > >;

  double const initialPrimarySpeciesConcentration[numPrimarySpecies] =
  {
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89000, // Cl-
    1.65e-2, // Mg+2
    1.09000 // Na+1
  };



  double const logInitialPrimarySpeciesConcentration[numPrimarySpecies] =
  {
    log( initialPrimarySpeciesConcentration[0] ),
    log( initialPrimarySpeciesConcentration[1] ),
    log( initialPrimarySpeciesConcentration[2] ),
    log( initialPrimarySpeciesConcentration[3] ),
    log( initialPrimarySpeciesConcentration[4] ),
    log( initialPrimarySpeciesConcentration[5] ),
    log( initialPrimarySpeciesConcentration[6] )
  };

  double logPrimarySpeciesConcentration[numPrimarySpecies];
  EquilibriumReactionsType::enforceEquilibrium_LogAggregate( 0,
                                                             hpcReact::geochemistry::carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                             hpcReact::geochemistry::carbonateIdentityActivityParams,
                                                             logInitialPrimarySpeciesConcentration,
                                                             logPrimarySpeciesConcentration );

  double const expectedPrimarySpeciesConcentrations[numPrimarySpecies] =
  {
    0.00057358597611063442, // H+
    0.00029604457466591314, // HCO3-
    0.015486690880470029, // Ca+2
    0.0025362432534459978, // SO4-2
    1.8598530940823459, // Cl-
    0.0099729692401427945, // Mg+2
    1.0720360981956494 // Na+1
  };

  for( int r=0; r<numPrimarySpecies; ++r )
  {
    EXPECT_NEAR( exp( logPrimarySpeciesConcentration[r] ), expectedPrimarySpeciesConcentrations[r], 1.0e-8 * expectedPrimarySpeciesConcentrations[r] );
  }


}

//******************************************************************************
// The B-dot activity model, verified against EQ3NR.
//
// The nine complexation reactions solved here are the model EQ3NR solves, so the expected values
// are its converged molalities, read from eq36Database/carbonate.3o. See the README there.
//
// The parameters are carbonateNosolidActivityParamsEQ36 rather than carbonateNosolidActivityParams:
// the comparison is only meaningful with the EQ3/6 B-dot parameters, since the phreeqc set leaves
// several of these complexes without an ion size and so at gamma = 1.
//
// The tolerance is set by the reference, not by the solver. EQ3NR prints five significant figures,
// and it fits the Debye-Huckel A and B over its temperature grid rather than reading the 25 C
// entry, so its effective values differ from the patched ones in the fifth decimal. Every primary
// species here agrees to 8.9e-5 or better.
TEST( testEquilibriumReactions, testcarbonateSystem_Bdot )
{
  using namespace hpcReact::geochemistry;

  static constexpr int numPrimarySpecies = carbonateSystemType::numPrimarySpecies();

  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< double, int, int,
                                                                           carbonateNosolidActivityType >;

  double const aggregatePrimarySpeciesConcentration[numPrimarySpecies] =
  {
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89, // Cl-
    1.65e-2, // Mg+2
    1.09 // Na+1
  };

  // The totals themselves, which is what a caller has before a run. B-dot cannot be started from
  // them directly; enforceEquilibrium_Aggregate seeds itself with an ideal solve to get there.
  double logInitialGuess[numPrimarySpecies];
  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    logInitialGuess[i] = log( aggregatePrimarySpeciesConcentration[i] );
  }

  double logPrimarySpeciesConcentration[numPrimarySpecies];
  EquilibriumReactionsType::enforceEquilibrium_Aggregate( 298.15,
                                                          carbonateSystem.equilibriumReactionsParameters(),
                                                          carbonateNosolidActivityParamsEQ36,
                                                          aggregatePrimarySpeciesConcentration,
                                                          logInitialGuess,
                                                          logPrimarySpeciesConcentration );

  // EQ3NR converged molalities. This solve reports only the primary species; the secondary species
  // of the same run are compared in testCarbonateActivityVsEQ36.
  double const expectedPrimarySpeciesConcentrations[numPrimarySpecies] =
  {
    6.5867e-04, // H+
    6.1144e-04, // HCO3-
    3.2573e-02, // Ca+2
    1.4996e-02, // SO4-2
    1.8836e+00, // Cl-
    1.4435e-02, // Mg+2
    1.0766e+00 // Na+1
  };

  double const eq36Tolerance = 5.0e-4;

  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    EXPECT_NEAR( exp( logPrimarySpeciesConcentration[i] ),
                 expectedPrimarySpeciesConcentrations[i],
                 eq36Tolerance * expectedPrimarySpeciesConcentrations[i] );
  }
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
