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
TEST( testEquilibriumReactions, testcarbonateSystemAllEquilibrium )
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
  { 2.1579694253441686e-11, // OH-
    0.3755789961165058, // CO2
    3.4214835005538611e-11, // CO3-2
    1.9160014879413049e-05, // CaHCO3+
    0.0025013721967110923, // CaSO4
    0.030083903919781853, // CaCl+
    0.0028032598559295028, // CaCl2
    0.0063383337566393907, // MgSO4
    0.019567697287351467, // NaSO4-
    4.754873524374959e-05, // CaCO3
    0.00046855267453226469, // H+
    0.00035429509915661095, // HCO3-
    0.0032447552774548935, // Ca+2
    0.0036925967592983458, // SO4-2
    1.8543095763683592, // Cl-
    0.01016166624336071, // Mg+2
    1.0704323027126488 // Na+1
  };

  std::cout<<" RESIDUAL_FORM 0:"<<std::endl;
  testEnforceEquilibrium< double, 0 >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                       initialSpeciesConcentration,
                                       expectedSpeciesConcentrations );

  // std::cout<<" RESIDUAL_FORM 1:"<<std::endl;
  // testEnforceEquilibrium< double, 1 >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
  //                                      initialSpeciesConcentration,
  //                                      expectedSpeciesConcentrations );

  std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
  testEnforceEquilibrium< double, 2 >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                       initialSpeciesConcentration,
                                       expectedSpeciesConcentrations );

}


TEST( testEquilibriumReactions, testcarbonateSystemAllEquilibrium2 )
{

  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< double,
                                                                           int,
                                                                           int >;

  constexpr int numPrimarySpecies = hpcReact::geochemistry::carbonateSystemAllEquilibrium.numPrimarySpecies();

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
                                                             logInitialPrimarySpeciesConcentration,
                                                             logPrimarySpeciesConcentration );

  double const expectedPrimarySpeciesConcentrations[numPrimarySpecies] =
  {
    0.00046855267453254149, // H+
    0.00035429509915645743, // HCO3-
    0.0032447552774548518, // Ca+2
    0.0036925967592983211, // SO4-2
    1.8543095763683592, // Cl-
    0.010161666243360675, // Mg+2
    1.0704323027126488 // Na+1
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
