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

  // brine source
  double const initialAggregateSpeciesConcentration[numPrimarySpecies] =
  {
    5.50E-06, // H+
    2.75E-01, // Ca+2
    1.21E-01, // Mg+2
    3.41E+00, // Na+
    3.25E-02, // K+
    5.45E-04, // Al+++
    2.60E+00, // HCO3-
    3.81E-01, // SO4-2
    1.00E+00, // Cl-
    1.03E-02 //SiO2(aq)
  };
  //brine, from sink at t=0.23h
  //   double const initialAggregateSpeciesConcentration[numPrimarySpecies] =
  // {
  //   0.01202, // H+
  //   9.31E-02, // Ca+2
  //   1.36E-01, // Mg+2
  //   3.40, // Na+
  //   2.38E-01, // K+
  //   2.93E-03, // Al+++
  //   2.60, // HCO3-
  //   3.80E-01, // SO4-2
  //   1.00E+00, // Cl-
  //   2.79E-01 //SiO2(aq)
  // };

  double const expectedSpeciesConcentrations[numPrimarySpecies] =
  {
    1.9964e-04, // H+
    0.1050, // Ca+2
    0.1209, // Mg+2
    3.4000, // Na+
    0.1250, // K+
    5.4367e-04, // Al+++
    2.60, // HCO3-
    0.380, // SO4-2
    1.0, // Cl-
    0.0103 //SiO2(aq)
  };

  timeStepTest< double, true >( forgeSystem,
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
