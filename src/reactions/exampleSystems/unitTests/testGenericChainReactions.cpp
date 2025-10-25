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

#include "reactions/unitTestUtilities/kineticReactionsTestUtilities.hpp"
#include "../ChainGeneric.hpp"


using namespace hpcReact;
using namespace hpcReact::unitTest_utilities;

//******************************************************************************
TEST( testChainGenericKineticReactions, computeReactionRatesTest_chainReactionParams )
{
  using namespace hpcReact::ChainGeneric;

  double const initialSpeciesConcentration[] =
  {
    1.0, // C1
    1e-8, // C2
    1e-8 // C3
  };
  double const surfaceArea[] =
  {
    0.0, // C1 -> C2
    0.0, // C2 -> C3
    0.0 // C3 ->
  };

  double const expectedReactionRates[] =
  {
    0.05, // C1 -> C2
    3e-10, // C2 -> C3
    2e-10 // C3 ->
  };
  double const expectedReactionRatesDerivatives[][3] =
  { { 0.05, 0.0, 0.0 },
    { 0.0, 0.03, 0.0 },
    { 0.0, 0.0, 0.02 } };
  computeReactionRatesTest< double, false >( serialAllKineticParams.kineticReactionsParameters(),
                                             initialSpeciesConcentration,
                                             surfaceArea, // No use. Just to pass something here
                                             expectedReactionRates,
                                             expectedReactionRatesDerivatives );
  computeReactionRatesTest< double, true >( serialAllKineticParams.kineticReactionsParameters(),
                                            initialSpeciesConcentration,
                                            surfaceArea, // No use. Just to pass something here
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
