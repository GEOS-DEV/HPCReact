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
    0.00040311656239679382, // H+
    0.00041180885982392148, // HCO3-
    0.0032499045666604504, // Ca+2
    0.0036920967945592146, // SO4-2
    1.8542541730074311, // Cl-
    0.010162194793470079, // Mg+2
    1.070434904554991 // Na+1
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
