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

#include "constitutive/activity/Bdot.hpp"
#include "constitutive/activity/Identity.hpp"
#include "constitutive/ionicStrength/SpeciatedIonicStrength.hpp"

#include <gtest/gtest.h>


using namespace hpcReact;
using namespace hpcReact::reactionsSystems;
using namespace hpcReact::unitTest_utilities;

//******************************************************************************
TEST( testKineticReactions, computeReactionRatesTest_simpleKineticTestRateParams_Identity )
{
  using IonicStrengthType = bulkGeneric::simpleIonicStrengthType;
  using ActivityType = Identity< double, int, IonicStrengthType >;

  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const surfaceArea[] = { 0.0, 0.0 };

  double const expectedReactionRates[] = { 1.0, 0.25 };
  double const expectedReactionRatesDerivatives[][5] =
  { { 2.0, -0.5, 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 0.5, 0.25, 0.0 } };

  computeReactionRatesTest< double,
                            false,
                            ActivityType >( bulkGeneric::simpleKineticTestRateParams,
                                            bulkGeneric::simpleIdentityActivityTestParams,
                                            initialSpeciesConcentration,
                                            surfaceArea, // No use. Just to pass something here
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
  computeReactionRatesTest< double,
                            true,
                            ActivityType >( bulkGeneric::simpleKineticTestRateParams,
                                            bulkGeneric::simpleIdentityActivityTestParams,
                                            initialSpeciesConcentration,
                                            surfaceArea, // No use. Just to pass something here
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
}


TEST( testKineticReactions, computeReactionRatesTest_simpleKineticTestRateParams_Bdot )
{
  using IonicStrengthType = bulkGeneric::simpleIonicStrengthType;
  using ActivityType = Bdot< double, int, IonicStrengthType >;

  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const surfaceArea[] = { 0.0, 0.0 };

  // Reference values derived independently from the B-dot equations in
  // bdotKineticReference.py (I = 2.25, gamma = { 0.157531, 0.880784, 0.880784, 1,
  // 0.880784 }) and cross-checked against a central difference of the rates.
  double const expectedReactionRates[] = { 0.024816095046886668, 0.22019609961589134 };
  double const expectedReactionRatesDerivatives[][5] =
  { { 0.054908038968070227, -0.38657161606983814, 0.0013189622185741268, 0.0, 0.0013189622185742044 },
    { 0.078220259582388221, 0.019555064895597055, 0.45994726412737974, 0.22019609961589134, 0.019555064895596979 } };

  computeReactionRatesTest< double,
                            false,
                            ActivityType >( bulkGeneric::simpleKineticTestRateParams,
                                            bulkGeneric::simpleActivityTestParams,
                                            initialSpeciesConcentration,
                                            surfaceArea, // No use. Just to pass something here
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
  computeReactionRatesTest< double,
                            true,
                            ActivityType >( bulkGeneric::simpleKineticTestRateParams,
                                            bulkGeneric::simpleActivityTestParams,
                                            initialSpeciesConcentration,
                                            surfaceArea, // No use. Just to pass something here
                                            expectedReactionRates,
                                            expectedReactionRatesDerivatives );
}


TEST( testKineticReactions, computeSpeciesRatesTest_simpleKineticTestRateParams_Identity )
{
  using IonicStrengthType = bulkGeneric::simpleIonicStrengthType;
  using ActivityType = Identity< double, int, IonicStrengthType >;

  double const initialSpeciesConcentration[5] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedSpeciesRates[5] = { -2.0, 1.0, 0.75, -0.25, 0.5 };
  double const expectedSpeciesRatesDerivatives[5][5] = { { -4.0, 1.0, 0.0, 0.0, 0.0 },
    {  2.0, -0.5, 0.0, 0.0, 0.0 },
    {  2.0, -0.5, -0.5, -0.25, 0.0 },
    {  0.0, 0.0, -0.5, -0.25, 0.0 },
    {  0.0, 0.0, 1.0, 0.5, 0.0 } };

  computeSpeciesRatesTest< double,
                           false,
                           ActivityType >( bulkGeneric::simpleKineticTestRateParams,
                                           bulkGeneric::simpleIdentityActivityTestParams,
                                           initialSpeciesConcentration,
                                           expectedSpeciesRates,
                                           expectedSpeciesRatesDerivatives );

  computeSpeciesRatesTest< double,
                           true,
                           ActivityType >( bulkGeneric::simpleKineticTestRateParams,
                                           bulkGeneric::simpleIdentityActivityTestParams,
                                           initialSpeciesConcentration,
                                           expectedSpeciesRates,
                                           expectedSpeciesRatesDerivatives );

}


// TEST( testKineticReactions, testTimeStep )
// {
//   double const initialSpeciesConcentration[5] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
//   double const expectedSpeciesConcentrations[5] = { 3.92138293924124e-01, 3.03930853037938e-01, 5.05945480771998e-01,
// 7.02014627734060e-01, 5.95970744531880e-01 };

//   timeStepTest< double, false >( bulkGeneric::simpleKineticTestRateParams,
//                                  2.0,
//                                  10,
//                                  initialSpeciesConcentration,
//                                  expectedSpeciesConcentrations );

// }

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
