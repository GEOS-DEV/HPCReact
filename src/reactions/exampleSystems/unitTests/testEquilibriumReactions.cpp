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
#include "../BulkGeneric.hpp"

#include <gtest/gtest.h>

using namespace hpcReact;
using namespace hpcReact::unitTest_utilities;

//******************************************************************************
TEST( testEquilibriumReactions, computeResidualAndJacobianTest )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };


  {
    std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
    double const expectedResiduals[] = { -37.534508668465, -72.989575795250 };
    double const expectedJacobian[2][2] =
    { { 1.0e16, -2.0 },
      { -2.0, 4.0e16 } };
    computeResidualAndJacobianTest< double, 2 >( bulkGeneric::simpleTestRateParams,
                                                 initialSpeciesConcentration,
                                                 expectedResiduals,
                                                 expectedJacobian );
  }

}

//******************************************************************************
TEST( testEquilibriumReactions, testEnforceEquilibrium )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedSpeciesConcentrations[5] = { 3.92138294e-01, 3.03930853e-01, 5.05945481e-01, 7.02014628e-01, 5.95970745e-01 };


  std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
  testEnforceEquilibrium< double, 2 >( bulkGeneric::simpleTestRateParams.equilibriumReactionsParameters(),
                                       initialSpeciesConcentration,
                                       expectedSpeciesConcentrations );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
