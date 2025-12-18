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

#if defined(__INTELLISENSE__)
#include "EquilibriumReactions.hpp"
#endif

#include "reactions/massActions/MassActions.hpp"

namespace hpcReact
{
namespace reactionsSystems
{


template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE >::computeResidualAndJacobianAggregatePrimaryConcentrations( RealType const & temperature,
                                                                                              PARAMS_DATA const & params,
                                                                                              ARRAY_1D_TO_CONST const & targetAggregatePrimaryConcentrations,
                                                                                              ARRAY_1D_TO_CONST2 const & logPrimarySpeciesConcentration,
                                                                                              ARRAY_1D & residual,
                                                                                              ARRAY_2D & jacobian )
{
  HPCREACT_UNUSED_VAR( temperature );
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  RealType aggregatePrimaryConcentrations[numPrimarySpecies] = {0.0};
  ARRAY_2D dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations = {{{0.0}}};
  massActions::calculateAggregatePrimaryConcentrationsWrtLogC< REAL_TYPE, INT_TYPE, INDEX_TYPE >( params,
                                                                                                  logPrimarySpeciesConcentration,
                                                                                                  aggregatePrimaryConcentrations,
                                                                                                  dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );


  for( IndexType i=0; i<numPrimarySpecies; ++i )
  {
    residual[i] = -(1.0 - aggregatePrimaryConcentrations[i] / targetAggregatePrimaryConcentrations[i]);
    for( IndexType j=0; j<numPrimarySpecies; ++j )
    {
      jacobian( i, j ) = -dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations[i][j] / targetAggregatePrimaryConcentrations[i];
    }
  }
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST >
HPCREACT_HOST_DEVICE inline
void
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE >::enforceEquilibrium_LogAggregate( REAL_TYPE const & temperature,
                                                                     PARAMS_DATA const & params,
                                                                     ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration0,
                                                                     ARRAY_1D & logPrimarySpeciesConcentration )
{
  HPCREACT_UNUSED_VAR( temperature );
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();
  double targetAggregatePrimarySpeciesConcentration[numPrimarySpecies] = { 0.0 };



  for( int i=0; i<numPrimarySpecies; ++i )
  {
    targetAggregatePrimarySpeciesConcentration[i] = logmath::exp( logPrimarySpeciesConcentration0[i] );
  }

  enforceEquilibrium_Aggregate( temperature,
                                params,
                                targetAggregatePrimarySpeciesConcentration,
                                logPrimarySpeciesConcentration0,
                                logPrimarySpeciesConcentration );

}


template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST >
HPCREACT_HOST_DEVICE inline
void
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE >::enforceEquilibrium_Aggregate( REAL_TYPE const & temperature,
                                                                  PARAMS_DATA const & params,
                                                                  ARRAY_1D_TO_CONST const & targetAggregatePrimarySpeciesConcentration,
                                                                  ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration0,
                                                                  ARRAY_1D & logPrimarySpeciesConcentration )
{
  if constexpr ( PARAMS_DATA::numSecondarySpecies() <= 0 )
  {
    return;
  }

  HPCREACT_UNUSED_VAR( temperature );
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  double residual[numPrimarySpecies] = { 0.0 };
//  double aggregatePrimarySpeciesConcentration[numPrimarySpecies] = { 0.0 };
  double dLogCp[numPrimarySpecies] = { 0.0 };
  CArrayWrapper< double, numPrimarySpecies, numPrimarySpecies > jacobian;

  for( int i=0; i<numPrimarySpecies; ++i )
  {
    logPrimarySpeciesConcentration[i] = logPrimarySpeciesConcentration0[i];
  }


  REAL_TYPE residualNorm = 0.0;
  // // Print for MoMaS only
  // //         0:     1e-20       -0           2 -2.5e+11       1e-20        7           2      1.8           1        5
  // printf( "iter       X1       R0           X2      R1          X3       R2          X4       R3           S       R4\n" );
  // printf( "----   ---------------      ---------------      ---------------      ---------------      ---------------\n" );
  for( int k=0; k<150; ++k )
  {
    computeResidualAndJacobianAggregatePrimaryConcentrations( temperature,
                                                              params,
                                                              targetAggregatePrimarySpeciesConcentration,
                                                              logPrimarySpeciesConcentration,
                                                              residual,
                                                              jacobian );

    residualNorm = 0.0;
    for( int i = 0; i < numPrimarySpecies; ++i )
    {
      residualNorm += residual[i] * residual[i];
    }
    residualNorm = sqrt( residualNorm );

    //  // Print for MoMaS only
    // printf( "%2d:  %8.2g %8.2g    %8.2g %8.2g    %8.2g %8.2g    %8.2g %8.2g    %8.2g %8.2g  \n",
    //         k,
    //         logmath::exp( logPrimarySpeciesConcentration[0] ),
    //         residual[0],
    //         logmath::exp( logPrimarySpeciesConcentration[1] ),
    //         residual[1],
    //         logmath::exp( logPrimarySpeciesConcentration[2] ),
    //         residual[2],
    //         logmath::exp( logPrimarySpeciesConcentration[3] ),
    //         residual[3],
    //         logmath::exp( logPrimarySpeciesConcentration[4] ),
    //         residual[4] );

    //printf( "iter, residualNorm = %2d, %16.10g \n", k, residualNorm );
    if( residualNorm < 1.0e-12 )
    {
      printf( " converged\n" );
      break;
    }

    solveNxN_pivoted< double, numPrimarySpecies >( jacobian.data, residual, dLogCp );


    for( IndexType i=0; i<numPrimarySpecies; ++i )
    {
      logPrimarySpeciesConcentration[i] += dLogCp[i];
    }

  }
}

} // namespace reactionsSystems
} // namespace hpcReact
