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

#include "constitutive/activity/Bdot.hpp"
#include "constitutive/activity/Identity.hpp"

#include <gtest/gtest.h>

using namespace hpcReact;
using namespace hpcReact::unitTest_utilities;

using SimpleBdotActivityType     = Bdot< double, int, bulkGeneric::simpleIonicStrengthType >;
using SimpleIdentityActivityType = Identity< double, int, bulkGeneric::simpleIonicStrengthType >;

//******************************************************************************
TEST( testEquilibriumReactions, computeResidualAndJacobianTest_Identity )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };


  {
    std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
    double const expectedResiduals[] = { -37.534508668465, -72.989575795250 };
    double const expectedJacobian[2][2] =
    { { 1.0e16, -2.0 },
      { -2.0, 4.0e16 } };

    computeResidualAndJacobianTest< double, 2, SimpleIdentityActivityType >( bulkGeneric::simpleTestRateParams,
                                                                             bulkGeneric::simpleIdentityActivityTestParams,
                                                                             initialSpeciesConcentration,
                                                                             expectedResiduals,
                                                                             expectedJacobian );
  }

}

//******************************************************************************
TEST( testEquilibriumReactions, computeResidualAndJacobianTest_Bdot )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };


  {
    std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
    // Reference values from bdotEquilibriumExtentReference.py
    double const expectedResiduals[] = { -34.092130639267268, -73.116518201896085 };
    double const expectedJacobian[2][2] =
    { { 1.0e16, -1.8755344687582751 },
      { -2.5328449939769717, 4.0e16 } };

    computeResidualAndJacobianTest< double, 2, SimpleBdotActivityType >( bulkGeneric::simpleTestRateParams,
                                                                         bulkGeneric::simpleActivityTestParams,
                                                                         initialSpeciesConcentration,
                                                                         expectedResiduals,
                                                                         expectedJacobian );
  }

}

//******************************************************************************
TEST( testEquilibriumReactions, testEnforceEquilibrium_Identity )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedSpeciesConcentrations[5] = { 3.92138294e-01, 3.03930853e-01, 5.05945481e-01, 7.02014628e-01, 5.95970745e-01 };


  std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
  testEnforceEquilibrium< double, 2, SimpleIdentityActivityType >( bulkGeneric::simpleTestRateParams.equilibriumReactionsParameters(),
                                                                   bulkGeneric::simpleIdentityActivityTestParams,
                                                                   initialSpeciesConcentration,
                                                                   expectedSpeciesConcentrations );

}

//******************************************************************************
TEST( testEquilibriumReactions, testEnforceEquilibrium_Bdot )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  // Reference values from bdotEquilibriumExtentReference.py.
  double const expectedSpeciesConcentrations[5] =
  { 0.84988563454963306, 0.075057182725183552, 0.31542526797552606, 0.74036808525034259, 0.51926382949931493 };


  std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
  testEnforceEquilibrium< double, 2, SimpleBdotActivityType >( bulkGeneric::simpleTestRateParams.equilibriumReactionsParameters(),
                                                               bulkGeneric::simpleActivityTestParams,
                                                               initialSpeciesConcentration,
                                                               expectedSpeciesConcentrations );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
