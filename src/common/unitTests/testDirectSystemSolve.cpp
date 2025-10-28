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


#include "../DirectSystemSolve.hpp"
#include "common/pmpl.hpp"

#include <gtest/gtest.h>

using namespace hpcReact;

template< typename REAL_TYPE, int N >
struct LinearSystem
{
  REAL_TYPE A[N][N];
  REAL_TYPE b[N];
  REAL_TYPE x[N];
};


void test3x3_helper()
{
  // **Define a Sample NxN Linear System**

  LinearSystem< double, 3 > linearSystem
  {
    { { 1.0, 2.0, 3.0 },
      { 2.0, -1.0, 1.0 },
      { 3.0, 4.0, 5.0 }
    },
    { 14.0, 3.0, 24.0 }, // Right-hand side
    { 0.0, 0.0, 0.0 } // Solution
  };

  pmpl::genericKernelWrapper( 1, &linearSystem, [] HPCREACT_DEVICE ( auto * const copyOfLinearSystem )
  {
    solveNxN_pivoted< double, 3 >( copyOfLinearSystem->A, copyOfLinearSystem->b, copyOfLinearSystem->x );
  } );

  EXPECT_NEAR( linearSystem.x[0], 0.0, std::numeric_limits< double >::epsilon()*100 );
  EXPECT_NEAR( linearSystem.x[1], 1.0, std::numeric_limits< double >::epsilon()*100 );
  EXPECT_NEAR( linearSystem.x[2], 4.0, std::numeric_limits< double >::epsilon()*100 );

}

TEST( testDirectSystemSolve, test3x3 )
{
  test3x3_helper();
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
