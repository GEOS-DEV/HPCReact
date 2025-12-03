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


TEST( testForgeReactions, testTimeStep_forgeSystem )
{
  using namespace hpcReact::geochemistry;

  static constexpr int numPrimarySpecies = forgeSystemType::numPrimarySpecies();

  double const surfaceArea[forgeSystemType::numKineticReactions()] =
  {
    1.0, 1.0, 1.0, 1.0, 1.0 // SiO2(s), KAlMg3Si3O10(OH)2(s), CaAl2(SiO4)2(s), KAlSi3O8(s), Al2Si2O5(OH)4(s)
  };

  double const initialAggregateSpeciesConcentration[numPrimarySpecies] =
  {
    5e-08, // H+
    8.3991e-0, // Ca+2
    0.00014146, // Mg+2
    0.0034061, // Na+
    0.00020949, // K+
    2.3494e-08, // Al+++
    0.0016047, // HCO3-
    0.00038069, // SO4-2
    0.0010018, // Cl-
    4.4047e-06 //SiO2(aq)
  };

  double const expectedSpeciesConcentrations[numPrimarySpecies] =
  {
    5e-08, // H+
    8.3991e-0, // Ca+2
    0.00014146, // Mg+2
    0.0034061, // Na+
    0.00020949, // K+
    2.3494e-08, // Al+++
    0.0016047, // HCO3-
    0.00038069, // SO4-2
    0.0010018, // Cl-
    4.4047e-06 //SiO2(aq)
  };

  timeStepTest< double, true >( forgeSystem,
                                1.e-3,
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
