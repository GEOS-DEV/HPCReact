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
#include "../Momas.hpp"

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

TEST( testEquilibriumReactions, testMoMasAllEquilibrium )
{

  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< double,
                                                                           int,
                                                                           int >;

  constexpr int numPrimarySpecies = hpcReact::geochemistry::momasSystemAllEquilibrium.numPrimarySpecies();

  double const targetAggregatePrimarySpeciesConcentration[numPrimarySpecies] =
  {
    1.0e-20,     // X1
    -2.0,    // X2
    1.0e-20,     // X3
    2.0,     // X4
    1.0      // S
  };

  double const initialPrimarySpeciesConcentration[numPrimarySpecies] =
  {
    1.0e-20,  // X1
    0.02,     // X2
    1.0e-20,  // X3
    1.0,     // X4
    1.00      // S
  };

  double const logInitialPrimarySpeciesConcentration[numPrimarySpecies] =
  {
    log( initialPrimarySpeciesConcentration[0] ),
    log( initialPrimarySpeciesConcentration[1] ),
    log( initialPrimarySpeciesConcentration[2] ),
    log( initialPrimarySpeciesConcentration[3] ),
    log( initialPrimarySpeciesConcentration[4] )
  };

  double logPrimarySpeciesConcentration[numPrimarySpecies];
  EquilibriumReactionsType::enforceEquilibrium_Aggregate( 0,
                                                          hpcReact::geochemistry::momasSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                          targetAggregatePrimarySpeciesConcentration,
                                                          logInitialPrimarySpeciesConcentration,
                                                          logPrimarySpeciesConcentration );

  double const expectedPrimarySpeciesConcentrations[numPrimarySpecies] =
  {
    9.7051090442170804E-21, // X1
    5.0023298955833342E-12, // X2
    1.9327372426296357E-33, // X3
    7.3929274619958745E-12, // X4
    9.8708294125907346E-13  // S
  };

  for( int r=0; r<numPrimarySpecies; ++r )
  {
    EXPECT_NEAR( exp( logPrimarySpeciesConcentration[r] ), expectedPrimarySpeciesConcentrations[r], 1.0e-8 * expectedPrimarySpeciesConcentrations[r] );
  }


}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}