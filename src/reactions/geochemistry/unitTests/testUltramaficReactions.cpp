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


TEST( testUltramaficReactions, testTimeStep_ultramaficSystem )
{
  using namespace hpcReact::geochemistry;

  static constexpr int numPrimarySpecies = ultramaficSystemType::numPrimarySpecies();

  double const surfaceArea[ultramaficSystemType::numKineticReactions()] =
  {
   1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 // Forsterite, Enstatite, Chrysotile, Diopside, Calcite, Dolomite-ord, Magnesite, Talc, Tremolite, Brucite
  };
  
  // Injection solution from Oman's paper (H+  Cl- HCO3- Ca++ Mg++ Na+ SO4-- SiO2(aq))
  double const initialAggregateSpeciesConcentration[numPrimarySpecies] =
  {
    7.94328E-05,  // H+    
    4.32000E-03,  //Cl                                                           
    122.30,  // HCO3- 
    1.73000E-03, // Ca++                                                               
    0.55E-03,    // Mg++   
    2.23000E-03, // Na+
    0.60000E-03, // SO4--                                                                                  
    0.00400E-03  // SiO2(aq)
  };

  // Solution of equilibrated groundwater (Oman's paper) with serpentine rock in EQ3/6
  //  double const initialAggregateSpeciesConcentration[numPrimarySpecies] =
  //  {
  //   9.55e-12,  // H+    
  //   8.90262E-03,  //Cl                                                           
  //   9.15090E-06,  // HCO3- 
  //   5.85434E-03, // Ca++                                                               
  //   6.88654E-08,    // Mg++   
  //   1.55841E-02, // Na+
  //   5.83275E-04, // SO4--                                                                                  
  //   1.78572E-06  // SiO2(aq)
  //  };

  // Currently the expected concentrations are random values
  double const expectedSpeciesConcentrations[numPrimarySpecies] =
  {
    7.94328E-05,  // H+    
    4.32000E-03,  //Cl                                                           
    1.22300E-01,  // HCO3- 
    1.73000E-03, // Ca++                                                               
    0.55E-03,    // Mg++   
    2.23000E-03, // Na+
    0.60000E-03, // SO4--                                                                                  
    0.00400E-03  // SiO2(aq)
  };

  timeStepTest< double, true >( ultramaficSystem,
                                1.e-5,
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
