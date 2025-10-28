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


#include "reactions/reactionsSystems/EquilibriumReactions.hpp"
#include "common/macros.hpp"
#include "common/pmpl.hpp"

#include <gtest/gtest.h>


namespace hpcReact
{

namespace unitTest_utilities
{

template< typename REAL_TYPE >
REAL_TYPE tolerance( REAL_TYPE const a, REAL_TYPE const b )
{
  return std::numeric_limits< double >::epsilon() * std::max( fabs( a ), fabs( b ) ) * 10;
}

template< int numReactions, int numSpecies >
struct ComputeResidualAndJacobianTestData 
{
  CArrayWrapper< double, numReactions, numReactions > jacobian;
  double residual[numReactions] = { 0.0 };
  double speciesConcentration[numSpecies];
};

//******************************************************************************
template< typename REAL_TYPE,
          int RESIDUAL_FORM,
          typename PARAMS_DATA >
void computeResidualAndJacobianTest( PARAMS_DATA const & params,
                                     REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies()],
                                     REAL_TYPE const (&expectedResidual)[PARAMS_DATA::numReactions()],
                                     REAL_TYPE const (&expectedJacobian)[PARAMS_DATA::numReactions()][PARAMS_DATA::numReactions()] )
{
  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< REAL_TYPE,
                                                                           int,
                                                                           int >;

  static constexpr int numSpecies = PARAMS_DATA::numSpecies();
  static constexpr int numReactions = PARAMS_DATA::numReactions();

  double const temperature = 298.15;

  ComputeResidualAndJacobianTestData<numReactions, numSpecies> data;
  for( int i = 0; i < numSpecies; ++i )
  {
    data.speciesConcentration[i] = initialSpeciesConcentration[i];
  } 

  pmpl::genericKernelWrapper( 1, &data, [params, temperature] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    double xi[numReactions] = { 0.0 };

    EquilibriumReactionsType::computeResidualAndJacobianReactionExtents( temperature,
                                                                         params,
                                                                         dataCopy->speciesConcentration,
                                                                         xi,
                                                                         dataCopy->residual,
                                                                         dataCopy->jacobian );
  });

//  printf( "R = { %8.4g, %8.4g }\n", residual[0], residual[1] );
  for( int r=0; r<numReactions; ++r )
  {
    EXPECT_NEAR( data.residual[r], expectedResidual[r], 1.0e-8 );
  }

  // HPCREACT_UNUSED_VAR( expectedJacobian );
  // printf( "J = { { %8.4g, %8.4g },\n    { %8.4g %8.4g } }\n\n",
  //         jacobian( 0, 0 ), jacobian( 0, 1 ),
  //         jacobian( 1, 0 ), jacobian( 1, 1 ) );
  for( int a = 0; a < numReactions; ++a )
  {
    for( int b = 0; b < numReactions; ++b )
    {
      EXPECT_NEAR( data.jacobian( a, b ), expectedJacobian[a][b], tolerance( data.jacobian( a, b ), expectedJacobian[a][b] ) );
    }
  }
}


//******************************************************************************

template< int numSpecies >
struct TestEnforceEquilibriumData
{
  double speciesConcentration0[numSpecies];
  double speciesConcentration[numSpecies];
};

template< typename REAL_TYPE,
          int RESIDUAL_FORM,
          typename PARAMS_DATA >
void testEnforceEquilibrium( PARAMS_DATA const & params,
                             REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies()],
                             REAL_TYPE const (&expectedSpeciesConcentrations)[PARAMS_DATA::numSpecies()] )
{
  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< REAL_TYPE,
                                                                           int,
                                                                           int >;

  static constexpr int numSpecies = PARAMS_DATA::numSpecies();

  double const temperature = 298.15;

  TestEnforceEquilibriumData<numSpecies> data;
  for( int i = 0; i < numSpecies; ++i )
  {
    data.speciesConcentration0[i] = initialSpeciesConcentration[i];
  }

  pmpl::genericKernelWrapper( 1, &data, [params, temperature] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    EquilibriumReactionsType::enforceEquilibrium_Extents( temperature,
                                                          params,
                                                          dataCopy->speciesConcentration0,
                                                          dataCopy->speciesConcentration );
  });

  for( int r=0; r<numSpecies; ++r )
  {
//    printf( "c[%d] = %22.14e\n", r, speciesConcentration[r] );
    EXPECT_NEAR( data.speciesConcentration[r], expectedSpeciesConcentrations[r], 1.0e-8 * expectedSpeciesConcentrations[r] );
  }

}

}

} // namespace hpcReact
