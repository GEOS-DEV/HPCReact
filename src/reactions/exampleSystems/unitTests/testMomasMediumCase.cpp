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
#include "../MoMasBenchmark.hpp"

using namespace hpcReact;
using namespace hpcReact::MoMasBenchmark;
using namespace hpcReact::unitTest_utilities;

//******************************************************************************

TEST( testEquilibriumReactions, testMoMasMediumEquilibrium )
{
  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< double,
                                                                           int,
                                                                           int >;

  constexpr int numPrimarySpecies = hpcReact::MoMasBenchmark::mediumCaseParams.numPrimarySpecies();

  double const targetAggregatePrimarySpeciesConcentration[numPrimarySpecies] =
  {
    1.0e-20, // X1
    -3.0, // X2
    1.0e-20, // X3
    1.0, // X4
    1.0 // S
  };

  double const initialPrimarySpeciesConcentration[numPrimarySpecies] =
  {
    1.0e-20, // X1
    0.02, // X2
    1.0e-20, // X3
    1.0, // X4
    1.0 // S
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
                                                          hpcReact::MoMasBenchmark::mediumCaseParams.equilibriumReactionsParameters(),
                                                          targetAggregatePrimarySpeciesConcentration,
                                                          logInitialPrimarySpeciesConcentration,
                                                          logPrimarySpeciesConcentration );

  double const expectedPrimarySpeciesConcentrations[numPrimarySpecies] =
  {
    9.9999999999999919e-21, // X1
    0.14796989521717838, // X2
    5.7165444793692536e-24, // X3
    0.025616412699749774, // X4
    0.53958559521499294 // S
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
