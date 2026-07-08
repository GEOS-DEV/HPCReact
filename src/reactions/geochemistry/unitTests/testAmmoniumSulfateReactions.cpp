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


TEST( testMascagniteReactions, testTimeStep_mascagniteSystem )
{
  using namespace hpcReact::geochemistry;

  static constexpr int numPrimarySpecies = mascagniteSystemType::numPrimarySpecies();

  double const surfaceArea[mascagniteSystemType::numKineticReactions()] =
  {
    1.0 // (NH4)2SO4(s)
  };

  // Initial concentrations
  double const initialAggregateSpeciesConcentration[numPrimarySpecies] =
  {
    3.16e-6, // H+
    1.426147, // NH4+  
    1.0, // SO4-- 
  };


  double const expectedSpeciesConcentrations[numPrimarySpecies] =
  {
    1.9964e-04, // H+
    0.1050, // NH4+
    0.1209, // SO4--
  };

  timeStepTest< double, true >( mascagniteSystem,
                                1.e-7,
                                30,
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
