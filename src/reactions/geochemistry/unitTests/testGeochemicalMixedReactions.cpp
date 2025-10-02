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

#include "reactions/unitTestUtilities/mixedReactionsTestUtilities.hpp"
#include "../GeochemicalSystems.hpp"


using namespace hpcReact;
using namespace hpcReact::unitTest_utilities;


TEST( testMixedReactions, testTimeStep_carbonateSystem )
{
  using namespace hpcReact::geochemistry;

  constexpr int numPrimarySpecies = carbonateSystemType::numPrimarySpecies();

  double const surfaceArea[carbonateSystemType::numKineticReactions()] =
  {
    1.0, // CaCO3
  };

  double const initialAggregateSpeciesConcentration[numPrimarySpecies] =
  {
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89, // Cl-
    1.65e-2, // Mg+2
    1.09 // Na+1
  };

  double const expectedSpeciesConcentrations[numPrimarySpecies] =
  {
    0.00041832456577690664, // H+
    0.00039683803953351703, // HCO3-
    0.0032486426324025838, // Ca+2
    0.0036922193074066475, // SO4-2
    1.8542677503336003, // Cl-
    0.010162065270892553, // Mg+2
    1.0704342669906313 // Na+1
  };

  timeStepTest< double, true >( carbonateSystem,
                                1.0,
                                10,
                                initialAggregateSpeciesConcentration,
                                surfaceArea,
                                expectedSpeciesConcentrations );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
