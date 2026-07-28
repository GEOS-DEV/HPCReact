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

#pragma once

#include "common/macros.hpp"
#include <math.h>
#include <functional>
#include <iostream>

namespace hpcReact
{
namespace massActions
{

namespace massActions_impl
{

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename FUNC >
HPCREACT_HOST_DEVICE
inline
void calculateLogSecondaryActivities( PARAMS_DATA const & params,
                                      ARRAY_1D_TO_CONST const & logPrimaryActivities,
                                      ARRAY_1D & logSecondaryActivities,
                                      FUNC && derivativeFunc )
{
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  static constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();

  for( INDEX_TYPE i = 0; i < numSecondarySpecies; ++i )
  {
    logSecondaryActivities[i] = 0.0;
  }

  for( int j=0; j<numSecondarySpecies; ++j )
  {
    logSecondaryActivities[j] = -log( params.equilibriumConstant( j ) );
    for( int k=0; k<numPrimarySpecies; ++k )
    {
      logSecondaryActivities[j] += params.stoichiometricMatrix( j, k+numSecondarySpecies ) * ( logPrimaryActivities[k] );
      derivativeFunc( j, k, params.stoichiometricMatrix( j, k+numSecondarySpecies ) );
    }
  }
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D_TO_CONST,
          typename ARRAY_1D_SECONDARY,
          typename ARRAY_1D_PRIMARY,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateAggregatePrimaryConcentrationsFromActivitiesWrtLogC( PARAMS_DATA const & params,
                                                                   ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                                   ARRAY_1D_TO_CONST2 const & logPrimaryActivities,
                                                                   ARRAY_2D_TO_CONST const & dLogPrimaryActivities_dLogPrimarySpeciesConcentrations,
                                                                   ARRAY_1D_SECONDARY & logSecondaryActivities,
                                                                   ARRAY_1D_PRIMARY & aggregatePrimarySpeciesConcentrations,
                                                                   ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();

  for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
  {
    for( INDEX_TYPE j = 0; j < numPrimarySpecies; ++j )
    {
      dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations[i][j] = 0.0;
    }
  }

  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    REAL_TYPE const primarySpeciesConcentration_i = exp( logPrimarySpeciesConcentrations[i] );
    aggregatePrimarySpeciesConcentrations[i] = primarySpeciesConcentration_i;
    dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, i ) = primarySpeciesConcentration_i;
  }

  if constexpr( numSecondarySpecies <= 0 )
  {
    return;
  }

  massActions_impl::calculateLogSecondaryActivities< REAL_TYPE,
                                                     INT_TYPE,
                                                     INDEX_TYPE >( params,
                                                                   logPrimaryActivities,
                                                                   logSecondaryActivities,
                                                                   []( INDEX_TYPE, INDEX_TYPE, REAL_TYPE ){} );

  REAL_TYPE dLogSecondaryActivities_dLogPrimaryActivities[numSecondarySpecies][numPrimarySpecies] = {{ 0.0 }};
  massActions_impl::calculateLogSecondaryActivities< REAL_TYPE,
                                                     INT_TYPE,
                                                     INDEX_TYPE >( params,
                                                                   logPrimaryActivities,
                                                                   logSecondaryActivities,
                                                                   [&]( const int j, const int k, REAL_TYPE const value )
  {
    dLogSecondaryActivities_dLogPrimaryActivities[j][k] = value;
  } );

  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    for( int j = 0; j < numSecondarySpecies; ++j )
    {
      REAL_TYPE const secondarySpeciesConcentration_j = exp( logSecondaryActivities[j] );
      REAL_TYPE const nu_ji = params.stoichiometricMatrix( j, i + numSecondarySpecies );

      aggregatePrimarySpeciesConcentrations[i] += nu_ji * secondarySpeciesConcentration_j;

      for( int k = 0; k < numPrimarySpecies; ++k )
      {
        REAL_TYPE dLogSecondaryActivity_j_dLogPrimarySpeciesConcentration_k = 0.0;
        for( int l = 0; l < numPrimarySpecies; ++l )
        {
          dLogSecondaryActivity_j_dLogPrimarySpeciesConcentration_k +=
            dLogSecondaryActivities_dLogPrimaryActivities[j][l] *
            dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[l][k];
        }

        REAL_TYPE const dSecondarySpeciesConcentration_j_dLogPrimarySpeciesConcentration_k =
          secondarySpeciesConcentration_j * dLogSecondaryActivity_j_dLogPrimarySpeciesConcentration_k;

        dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, k ) +=
          nu_ji * dSecondarySpeciesConcentration_j_dLogPrimarySpeciesConcentration_k;
      }
    }
  }
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D_TO_CONST,
          typename ARRAY_1D_SECONDARY,
          typename ARRAY_1D_PRIMARY,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateTotalAndMobileAggregatePrimaryConcentrationsFromActivitiesWrtLogC( PARAMS_DATA const & params,
                                                                                 ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                                                 ARRAY_1D_TO_CONST2 const & logPrimaryActivities,
                                                                                 ARRAY_2D_TO_CONST const & dLogPrimaryActivities_dLogPrimarySpeciesConcentrations,
                                                                                 ARRAY_1D_SECONDARY & logSecondaryActivities,
                                                                                 ARRAY_1D_PRIMARY & aggregatePrimarySpeciesConcentrations,
                                                                                 ARRAY_1D_PRIMARY & mobileAggregatePrimarySpeciesConcentrations,
                                                                                 ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations,
                                                                                 ARRAY_2D & dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();

  for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
  {
    for( INDEX_TYPE j = 0; j < numPrimarySpecies; ++j )
    {
      dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations[i][j] = 0.0;
      dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations[i][j] = 0.0;
    }
  }

  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    REAL_TYPE const primarySpeciesConcentration_i = exp( logPrimarySpeciesConcentrations[i] );
    aggregatePrimarySpeciesConcentrations[i] = primarySpeciesConcentration_i;
    mobileAggregatePrimarySpeciesConcentrations[i] = primarySpeciesConcentration_i;
    dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, i ) = primarySpeciesConcentration_i;
    dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, i ) = primarySpeciesConcentration_i;
  }

  if constexpr( numSecondarySpecies <= 0 )
  {
    return;
  }

  massActions_impl::calculateLogSecondaryActivities< REAL_TYPE,
                                                     INT_TYPE,
                                                     INDEX_TYPE >( params,
                                                                   logPrimaryActivities,
                                                                   logSecondaryActivities,
                                                                   []( INDEX_TYPE, INDEX_TYPE, REAL_TYPE ){} );

  REAL_TYPE dLogSecondaryActivities_dLogPrimaryActivities[numSecondarySpecies][numPrimarySpecies] = {{ 0.0 }};
  massActions_impl::calculateLogSecondaryActivities< REAL_TYPE,
                                                     INT_TYPE,
                                                     INDEX_TYPE >( params,
                                                                   logPrimaryActivities,
                                                                   logSecondaryActivities,
                                                                   [&]( const int j, const int k, REAL_TYPE const value )
  {
    dLogSecondaryActivities_dLogPrimaryActivities[j][k] = value;
  } );

  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    for( int j = 0; j < numSecondarySpecies; ++j )
    {
      REAL_TYPE const secondarySpeciesConcentration_j = exp( logSecondaryActivities[j] );
      REAL_TYPE const nu_ji = params.stoichiometricMatrix( j, i + numSecondarySpecies );
      REAL_TYPE const mobileFlag = params.mobileSecondarySpeciesFlag( j );

      aggregatePrimarySpeciesConcentrations[i] += nu_ji * secondarySpeciesConcentration_j;
      mobileAggregatePrimarySpeciesConcentrations[i] += nu_ji * secondarySpeciesConcentration_j * mobileFlag;

      for( int k = 0; k < numPrimarySpecies; ++k )
      {
        REAL_TYPE dLogSecondaryActivity_j_dLogPrimarySpeciesConcentration_k = 0.0;
        for( int l = 0; l < numPrimarySpecies; ++l )
        {
          dLogSecondaryActivity_j_dLogPrimarySpeciesConcentration_k +=
            dLogSecondaryActivities_dLogPrimaryActivities[j][l] *
            dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[l][k];
        }

        REAL_TYPE const dSecondarySpeciesConcentration_j_dLogPrimarySpeciesConcentration_k =
          secondarySpeciesConcentration_j * dLogSecondaryActivity_j_dLogPrimarySpeciesConcentration_k;

        dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, k ) +=
          nu_ji * dSecondarySpeciesConcentration_j_dLogPrimarySpeciesConcentration_k;

        dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, k ) +=
          nu_ji * dSecondarySpeciesConcentration_j_dLogPrimarySpeciesConcentration_k * mobileFlag;
      }
    }
  }
}

} // namespace

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D >
HPCREACT_HOST_DEVICE
inline
void calculateLogSecondaryActivities( PARAMS_DATA const & params,
                                      ARRAY_1D_TO_CONST const & logPrimaryActivities,
                                      ARRAY_1D & logSecondaryActivities )
{
  if constexpr( PARAMS_DATA::numSecondarySpecies() <= 0 )
  {
    return;
  }

  massActions_impl::calculateLogSecondaryActivities< REAL_TYPE,
                                                     INT_TYPE,
                                                     INDEX_TYPE >( params,
                                                                   logPrimaryActivities,
                                                                   logSecondaryActivities,
                                                                   []( INDEX_TYPE, INDEX_TYPE, REAL_TYPE ){} );
}


template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateLogSecondaryActivitiesWrtLogC( PARAMS_DATA const & params,
                                             ARRAY_1D_TO_CONST const & logPrimaryActivities,
                                             ARRAY_1D & logSecondaryActivities,
                                             ARRAY_2D & dLogSecondaryActivities_dLogPrimaryActivities )
{
  massActions_impl::calculateLogSecondaryActivities< REAL_TYPE, INT_TYPE, INDEX_TYPE >( params,
                                                                                        logPrimaryActivities,
                                                                                        logSecondaryActivities,
                                                                                        [&]( const int j, const int k, REAL_TYPE const value )
  {
    dLogSecondaryActivities_dLogPrimaryActivities[j][k] = value;
  } );
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D >
HPCREACT_HOST_DEVICE
inline
void calculateLogSecondarySpeciesConcentration( PARAMS_DATA const & params,
                                                ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                ARRAY_1D & logSecondarySpeciesConcentrations )
{
  calculateLogSecondaryActivities< REAL_TYPE,
                                   INT_TYPE,
                                   INDEX_TYPE >( params,
                                                 logPrimarySpeciesConcentrations,
                                                 logSecondarySpeciesConcentrations );
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateLogSecondarySpeciesConcentrationWrtLogC( PARAMS_DATA const & params,
                                                       ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                       ARRAY_1D & logSecondarySpeciesConcentrations,
                                                       ARRAY_2D & dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations )
{
  calculateLogSecondaryActivitiesWrtLogC< REAL_TYPE,
                                          INT_TYPE,
                                          INDEX_TYPE >( params,
                                                        logPrimarySpeciesConcentrations,
                                                        logSecondarySpeciesConcentrations,
                                                        dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations );
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_PRIMARY,
          typename ARRAY_1D_SECONDARY,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateAggregatePrimaryConcentrationsWrtLogC( PARAMS_DATA const & params,
                                                     ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                     ARRAY_1D_SECONDARY & logSecondarySpeciesConcentrations,
                                                     ARRAY_1D_PRIMARY & aggregatePrimarySpeciesConcentrations,
                                                     ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  REAL_TYPE dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[numPrimarySpecies][numPrimarySpecies] = {{ 0.0 }};
  for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
  {
    dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[i][i] = 1.0;
  }

  massActions_impl::calculateAggregatePrimaryConcentrationsFromActivitiesWrtLogC< REAL_TYPE,
                                                                                  INT_TYPE,
                                                                                  INDEX_TYPE >( params,
                                                                                                logPrimarySpeciesConcentrations,
                                                                                                logPrimarySpeciesConcentrations,
                                                                                                dLogPrimaryActivities_dLogPrimarySpeciesConcentrations,
                                                                                                logSecondarySpeciesConcentrations,
                                                                                                aggregatePrimarySpeciesConcentrations,
                                                                                                dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateAggregatePrimaryConcentrationsWrtLogC( PARAMS_DATA const & params,
                                                     ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                     ARRAY_1D & aggregatePrimarySpeciesConcentrations,
                                                     ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();

  if constexpr( numSecondarySpecies > 0 )
  {
    REAL_TYPE logSecondarySpeciesConcentrations[numSecondarySpecies] = { 0.0 };
    calculateAggregatePrimaryConcentrationsWrtLogC< REAL_TYPE,
                                                    INT_TYPE,
                                                    INDEX_TYPE >( params,
                                                                  logPrimarySpeciesConcentrations,
                                                                  logSecondarySpeciesConcentrations,
                                                                  aggregatePrimarySpeciesConcentrations,
                                                                  dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
  }
  else
  {
    HPCREACT_UNUSED_VAR( logPrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( aggregatePrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
  }
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateAggregatePrimaryConcentrationsWrtLogC( PARAMS_DATA const & params,
                                                     ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                     ARRAY_1D_TO_CONST const & logPrimaryActivities,
                                                     ARRAY_1D & aggregatePrimarySpeciesConcentrations,
                                                     ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();

  if constexpr( numSecondarySpecies > 0 )
  {
    REAL_TYPE dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[numPrimarySpecies][numPrimarySpecies] = {{ 0.0 }};
    for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
    {
      dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[i][i] = 1.0;
    }

    REAL_TYPE logSecondarySpeciesConcentrations[numSecondarySpecies] = {0};

    massActions_impl::calculateAggregatePrimaryConcentrationsFromActivitiesWrtLogC< REAL_TYPE,
                                                                                    INT_TYPE,
                                                                                    INDEX_TYPE >( params,
                                                                                                  logPrimarySpeciesConcentrations,
                                                                                                  logPrimaryActivities,
                                                                                                  dLogPrimaryActivities_dLogPrimarySpeciesConcentrations,
                                                                                                  logSecondarySpeciesConcentrations,
                                                                                                  aggregatePrimarySpeciesConcentrations,
                                                                                                  dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
  }
  else
  {
    HPCREACT_UNUSED_VAR( logPrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( logPrimaryActivities );
    HPCREACT_UNUSED_VAR( aggregatePrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
  }
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateAggregatePrimaryConcentrationsWrtLogC( PARAMS_DATA const & params,
                                                     ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                     ARRAY_1D_TO_CONST2 const & logPrimaryActivities,
                                                     ARRAY_2D_TO_CONST const & dLogPrimaryActivities_dLogPrimarySpeciesConcentrations,
                                                     ARRAY_1D & aggregatePrimarySpeciesConcentrations,
                                                     ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();

  if constexpr( numSecondarySpecies > 0 )
  {
    REAL_TYPE logSecondarySpeciesConcentrations[numSecondarySpecies] = { 0.0 };
    massActions_impl::calculateAggregatePrimaryConcentrationsFromActivitiesWrtLogC< REAL_TYPE,
                                                                                    INT_TYPE,
                                                                                    INDEX_TYPE >( params,
                                                                                                  logPrimarySpeciesConcentrations,
                                                                                                  logPrimaryActivities,
                                                                                                  dLogPrimaryActivities_dLogPrimarySpeciesConcentrations,
                                                                                                  logSecondarySpeciesConcentrations,
                                                                                                  aggregatePrimarySpeciesConcentrations,
                                                                                                  dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
  }
  else
  {
    HPCREACT_UNUSED_VAR( logPrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( logPrimaryActivities );
    HPCREACT_UNUSED_VAR( dLogPrimaryActivities_dLogPrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( aggregatePrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
  }
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_PRIMARY,
          typename ARRAY_1D_SECONDARY,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateTotalAndMobileAggregatePrimaryConcentrationsWrtLogC( PARAMS_DATA const & params,
                                                                   ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                                   ARRAY_1D_SECONDARY & logSecondarySpeciesConcentrations,
                                                                   ARRAY_1D_PRIMARY & aggregatePrimarySpeciesConcentrations,
                                                                   ARRAY_1D_PRIMARY & mobileAggregatePrimarySpeciesConcentrations,
                                                                   ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations,
                                                                   ARRAY_2D & dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  REAL_TYPE dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[numPrimarySpecies][numPrimarySpecies] = {{ 0.0 }};
  for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
  {
    dLogPrimaryActivities_dLogPrimarySpeciesConcentrations[i][i] = 1.0;
  }

  massActions_impl::calculateTotalAndMobileAggregatePrimaryConcentrationsFromActivitiesWrtLogC< REAL_TYPE,
                                                                                                INT_TYPE,
                                                                                                INDEX_TYPE >( params,
                                                                                                              logPrimarySpeciesConcentrations,
                                                                                                              logPrimarySpeciesConcentrations,
                                                                                                              dLogPrimaryActivities_dLogPrimarySpeciesConcentrations,
                                                                                                              logSecondarySpeciesConcentrations,
                                                                                                              aggregatePrimarySpeciesConcentrations,
                                                                                                              mobileAggregatePrimarySpeciesConcentrations,
                                                                                                              dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations,
                                                                                                              dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D_TO_CONST,
          typename ARRAY_1D_PRIMARY,
          typename ARRAY_1D_SECONDARY,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateTotalAndMobileAggregatePrimaryConcentrationsWrtLogC( PARAMS_DATA const & params,
                                                                   ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                                   ARRAY_1D_TO_CONST2 const & logPrimaryActivities,
                                                                   ARRAY_2D_TO_CONST const & dLogPrimaryActivities_dLogPrimarySpeciesConcentrations,
                                                                   ARRAY_1D_SECONDARY & logSecondarySpeciesConcentrations,
                                                                   ARRAY_1D_PRIMARY & aggregatePrimarySpeciesConcentrations,
                                                                   ARRAY_1D_PRIMARY & mobileAggregatePrimarySpeciesConcentrations,
                                                                   ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations,
                                                                   ARRAY_2D & dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();

  if constexpr( numSecondarySpecies > 0 )
  {
    massActions_impl::calculateTotalAndMobileAggregatePrimaryConcentrationsFromActivitiesWrtLogC< REAL_TYPE,
                                                                                                  INT_TYPE,
                                                                                                  INDEX_TYPE >( params,
                                                                                                                logPrimarySpeciesConcentrations,
                                                                                                                logPrimaryActivities,
                                                                                                                dLogPrimaryActivities_dLogPrimarySpeciesConcentrations,
                                                                                                                logSecondarySpeciesConcentrations,
                                                                                                                aggregatePrimarySpeciesConcentrations,
                                                                                                                mobileAggregatePrimarySpeciesConcentrations,
                                                                                                                dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations,
                                                                                                                dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
  }
  else
  {
    HPCREACT_UNUSED_VAR( logPrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( logPrimaryActivities );
    HPCREACT_UNUSED_VAR( dLogPrimaryActivities_dLogPrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( logSecondarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( aggregatePrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( mobileAggregatePrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
  }
}

} // namespace massActions
} // namespace hpcReact
