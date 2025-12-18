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
#include "../BulkGeneric.hpp"

#include <gtest/gtest.h>


using namespace hpcReact;
using namespace hpcReact::reactionsSystems;
using namespace hpcReact::unitTest_utilities;

//******************************************************************************
TEST( testKineticReactions, computeReactionRatesTest_simpleKineticTestRateParams )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const surfaceArea[] = { 0.0, 0.0 };
  double const expectedReactionRates[] = { 1.0, 0.25 };
  double const expectedReactionRatesDerivatives[][5] =
  { { 2.0, -0.5, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.5, 0.25, 0.0 } };
  computeReactionRatesTest< double >( bulkGeneric::simpleKineticTestRateParams.kineticReactionsParameters(),
                                      initialSpeciesConcentration,
                                      surfaceArea,       // No use. Just to pass something here
                                      expectedReactionRates,
                                      expectedReactionRatesDerivatives );
}


TEST( testKineticReactions, computeSpeciesRatesTest_simpleKineticTestRateParams )
{
  double const initialSpeciesConcentration[5] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedSpeciesRates[5] = { -2.0, 1.0, 0.75, -0.25, 0.5 };
  double const expectedSpeciesRatesDerivatives[5][5] = { { -4.0, 1.0, 0.0, 0.0, 0.0 },
    {  2.0, -0.5, 0.0, 0.0, 0.0 },
    {  2.0, -0.5, -0.5, -0.25, 0.0 },
    {  0.0, 0.0, -0.5, -0.25, 0.0 },
    {  0.0, 0.0, 1.0, 0.5, 0.0 } };


  computeSpeciesRatesTest< double >( bulkGeneric::simpleKineticTestRateParams.kineticReactionsParameters(),
                                     initialSpeciesConcentration,
                                     expectedSpeciesRates,
                                     expectedSpeciesRatesDerivatives );

}


TEST( testKineticReactions, testTimeStep )
{
  double const initialSpeciesConcentration[5] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedSpeciesConcentrations[5] = { 3.92138293924124e-01, 3.03930853037938e-01, 5.05945480771998e-01, 7.02014627734060e-01, 5.95970744531880e-01 };

  timeStepTest< double >( bulkGeneric::simpleKineticTestRateParams.kineticReactionsParameters(),
                          2.0,
                          10,
                          initialSpeciesConcentration,
                          expectedSpeciesConcentrations );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
