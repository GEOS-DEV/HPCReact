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


#include "constitutive/activity/Bdot.hpp"
#include "reactions/reactionsSystems/Parameters.hpp"
#include "constitutive/ionicStrength/SpeciatedIonicStrength.hpp"
#include "common/pmpl.hpp"

#include <gtest/gtest.h>

using namespace hpcReact;




constexpr SpeciatedIonicStrength<double, int>::Params<3> testParams
{
  // Species charge
  { 1.0, -1.0, 2.0 }
};

TEST( testBdot, testIonicStrength )
{
  double speciesConcentration[ testParams.numSpecies() ] = { 0.1, 0.2, 0.3 };

  double I = SpeciatedIonicStrength< double, int >::calculate( testParams,
                                                               speciesConcentration );

  double expectedI = 0.5 * ( 0.1 * 1.0 * 1.0 + 0.2 * (-1.0) * (-1.0) + 0.3 * 2.0 * 2.0 );
  EXPECT_DOUBLE_EQ( I, expectedI );
  
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
