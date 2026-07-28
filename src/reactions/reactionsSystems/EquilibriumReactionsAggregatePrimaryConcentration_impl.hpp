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
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL >
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
                      INDEX_TYPE,
                      ACTIVITY_MODEL >::computeResidualAndJacobianAggregatePrimaryConcentrations( RealType const & temperature,
                                                                                                  PARAMS_DATA const & params,
                                                                                                  typename ACTIVITY_MODEL::Params const & activityParams,
                                                                                                  ARRAY_1D_TO_CONST const & targetAggregatePrimaryConcentrations,
                                                                                                  ARRAY_1D_TO_CONST2 const & logPrimarySpeciesConcentration,
                                                                                                  ARRAY_1D & residual,
                                                                                                  ARRAY_2D & jacobian )
{
  HPCREACT_UNUSED_VAR( temperature );
  static constexpr int numSpecies = PARAMS_DATA::numSpecies();
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  static constexpr int numSecondarySpeciesStorage = numSecondarySpecies > 0 ? numSecondarySpecies : 1;
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  RealType aggregatePrimaryConcentrations[numPrimarySpecies] = {0.0};
  RealType logPrimaryActivities[numPrimarySpecies] = {0.0};
  RealType dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[numPrimarySpecies][numPrimarySpecies] = {{0.0}};
  ARRAY_2D dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations = {{{0.0}}};

  if constexpr( numPrimarySpecies > 0 )
  {
    RealType logSpeciesConcentration[numSpecies] = {0.0};
    RealType logActivities[numSpecies] = {0.0};
    RealType dLogActivities_dLogSpeciesConcentrations[numSpecies][numSpecies] = {{0.0}};

    RealType dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[numSecondarySpeciesStorage][numPrimarySpecies] = {{0.0}};
    if constexpr( numSecondarySpecies > 0 )
    {
      RealType logSecondarySpeciesConcentrations[numSecondarySpecies] = {0.0};
      massActions::calculateLogSecondarySpeciesConcentrationWrtLogC< REAL_TYPE,
                                                                     INT_TYPE,
                                                                     INDEX_TYPE >( params,
                                                                                   logPrimarySpeciesConcentration,
                                                                                   logSecondarySpeciesConcentrations,
                                                                                   dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations );
      for( IndexType i = 0; i < numSecondarySpecies; ++i )
      {
        logSpeciesConcentration[i] = logSecondarySpeciesConcentrations[i];
      }
    }

    for( IndexType i = 0; i < numPrimarySpecies; ++i )
    {
      logSpeciesConcentration[i + numSecondarySpecies] = logPrimarySpeciesConcentration[i];
    }

    calculateActivities< RealType,
                         IntType,
                         IndexType,
                         ACTIVITY_MODEL,
                         true >( activityParams,
                                 logSpeciesConcentration,
                                 logActivities,
                                 dLogActivities_dLogSpeciesConcentrations );

    for( IndexType i = 0; i < numPrimarySpecies; ++i )
    {
      IndexType const fullRow = i + numSecondarySpecies;
      logPrimaryActivities[i] = logActivities[fullRow];

      for( IndexType j = 0; j < numPrimarySpecies; ++j )
      {
        RealType value = dLogActivities_dLogSpeciesConcentrations[fullRow][j + numSecondarySpecies];
        for( IndexType k = 0; k < numSecondarySpecies; ++k )
        {
          value += dLogActivities_dLogSpeciesConcentrations[fullRow][k] *
                   dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[k][j];
        }
        dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[i][j] = value;
      }
    }
  }

  massActions::calculateAggregatePrimaryConcentrationsWrtLogC< REAL_TYPE, INT_TYPE, INDEX_TYPE >( params,
                                                                                                  logPrimarySpeciesConcentration,
                                                                                                  logPrimaryActivities,
                                                                                                  dLogPrimaryActivities_dLogPrimarySpeciesConcentrations,
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
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST >
HPCREACT_HOST_DEVICE inline
void
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE,
                      ACTIVITY_MODEL >::enforceEquilibrium_LogAggregate( REAL_TYPE const & temperature,
                                                                         PARAMS_DATA const & params,
                                                                         typename ACTIVITY_MODEL::Params const & activityParams,
                                                                         ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration0,
                                                                         ARRAY_1D & logPrimarySpeciesConcentration )
{
  HPCREACT_UNUSED_VAR( temperature );
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();
  double targetAggregatePrimarySpeciesConcentration[numPrimarySpecies] = { 0.0 };



  for( int i=0; i<numPrimarySpecies; ++i )
  {
    targetAggregatePrimarySpeciesConcentration[i] = exp( logPrimarySpeciesConcentration0[i] );
  }

  enforceEquilibrium_Aggregate( temperature,
                                params,
                                activityParams,
                                targetAggregatePrimarySpeciesConcentration,
                                logPrimarySpeciesConcentration0,
                                logPrimarySpeciesConcentration );

}


template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST >
HPCREACT_HOST_DEVICE inline
void
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE,
                      ACTIVITY_MODEL >::enforceEquilibrium_Aggregate( REAL_TYPE const & temperature,
                                                                      PARAMS_DATA const & params,
                                                                      typename ACTIVITY_MODEL::Params const & activityParams,
                                                                      ARRAY_1D_TO_CONST const & targetAggregatePrimarySpeciesConcentration,
                                                                      ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration0,
                                                                      ARRAY_1D & logPrimarySpeciesConcentration )
{
  if constexpr( PARAMS_DATA::numSecondarySpecies() <= 0 )
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
                                                              activityParams,
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
    //         exp( logPrimarySpeciesConcentration[0] ),
    //         residual[0],
    //         exp( logPrimarySpeciesConcentration[1] ),
    //         residual[1],
    //         exp( logPrimarySpeciesConcentration[2] ),
    //         residual[2],
    //         exp( logPrimarySpeciesConcentration[3] ),
    //         residual[3],
    //         exp( logPrimarySpeciesConcentration[4] ),
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
