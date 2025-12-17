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
#include "common/logMath.hpp"
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
void calculateLogSecondarySpeciesConcentration( PARAMS_DATA const & params,
                                                ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                ARRAY_1D & logSecondarySpeciesConcentrations,
                                                FUNC && derivativeFunc )
{
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  static constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();

  for( INDEX_TYPE i = 0; i < numSecondarySpecies; ++i )
  {
    logSecondarySpeciesConcentrations[i] = 0.0;
  }

  for( int j=0; j<numSecondarySpecies; ++j )
  {
    REAL_TYPE const gamma = 1;
    logSecondarySpeciesConcentrations[j] = -logmath::log( params.equilibriumConstant( j ) ) - logmath::log( gamma );
    for( int k=0; k<numPrimarySpecies; ++k )
    {
      logSecondarySpeciesConcentrations[j] += params.stoichiometricMatrix( j, k+numSecondarySpecies ) * ( logPrimarySpeciesConcentrations[k] + logmath::log( gamma ) );
      derivativeFunc( j, k, params.stoichiometricMatrix( j, k+numSecondarySpecies ) );
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
void calculateLogSecondarySpeciesConcentration( PARAMS_DATA const & params,
                                                ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                ARRAY_1D & logSecondarySpeciesConcentrations )
{
  if constexpr( PARAMS_DATA::numSecondarySpecies() <= 0 )
  {
    return;
  }

  massActions_impl::calculateLogSecondarySpeciesConcentration< REAL_TYPE,
                                                               INT_TYPE,
                                                               INDEX_TYPE >( params,
                                                                             logPrimarySpeciesConcentrations,
                                                                             logSecondarySpeciesConcentrations,
                                                                             [](INDEX_TYPE, INDEX_TYPE, REAL_TYPE ){} );
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
  massActions_impl::calculateLogSecondarySpeciesConcentration< REAL_TYPE, INT_TYPE, INDEX_TYPE >( params,
                                                                                                  logPrimarySpeciesConcentrations,
                                                                                                  logSecondarySpeciesConcentrations,
                                                                                                  [&]( const int j, const int k, REAL_TYPE const value )
  {
    dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[j][k] = value;
  } );
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
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();


  calculateLogSecondarySpeciesConcentration< REAL_TYPE,
                                             INT_TYPE,
                                             INDEX_TYPE >( params,
                                                           logPrimarySpeciesConcentrations,
                                                           logSecondarySpeciesConcentrations );
  for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
  {
    for( INDEX_TYPE j = 0; j < numPrimarySpecies; ++j )
    {
      dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations[i][j] = 0.0;
    }
  }

  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    REAL_TYPE const speciesConcentration_i = logmath::exp( logPrimarySpeciesConcentrations[i] );
    aggregatePrimarySpeciesConcentrations[i] = speciesConcentration_i;
    dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, i ) = speciesConcentration_i;
    for( int j = 0; j < numSecondarySpecies; ++j )
    {
      REAL_TYPE const secondarySpeciesConcentrations_j = logmath::exp( logSecondarySpeciesConcentrations[j] );
      aggregatePrimarySpeciesConcentrations[i] += params.stoichiometricMatrix( j, i+numSecondarySpecies ) * secondarySpeciesConcentrations_j;
      for( int k=0; k<numPrimarySpecies; ++k )
      {
        REAL_TYPE const dSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentration = logmath::dWrtLogConst() * params.stoichiometricMatrix( j, k+numSecondarySpecies ) * secondarySpeciesConcentrations_j;
        dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, k ) += params.stoichiometricMatrix( j, i+numSecondarySpecies ) *
                                                                                                      dSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentration;
      }
    }
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
                                                     ARRAY_1D & aggregatePrimarySpeciesConcentrations,
                                                     ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();

  if constexpr( numSecondarySpecies > 0 )
  {
    REAL_TYPE logSecondarySpeciesConcentrations[numSecondarySpecies] = {0};

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
    GEOS_UNUSED_VAR( logPrimarySpeciesConcentrations, aggregatePrimarySpeciesConcentrations, dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
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
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();

  calculateLogSecondarySpeciesConcentration< REAL_TYPE,
                                             INT_TYPE,
                                             INDEX_TYPE >( params,
                                                           logPrimarySpeciesConcentrations,
                                                           logSecondarySpeciesConcentrations );
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
    REAL_TYPE const speciesConcentration_i = logmath::exp( logPrimarySpeciesConcentrations[i] );
    aggregatePrimarySpeciesConcentrations[i] = speciesConcentration_i;
    mobileAggregatePrimarySpeciesConcentrations[i] = speciesConcentration_i;
    dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, i ) = speciesConcentration_i;
    dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, i ) = speciesConcentration_i;

    for( int j = 0; j < numSecondarySpecies; ++j )
    {
      REAL_TYPE const secondarySpeciesConcentrations_j = logmath::exp( logSecondarySpeciesConcentrations[j] );
      aggregatePrimarySpeciesConcentrations[i] += params.stoichiometricMatrix( j, i+numSecondarySpecies ) * secondarySpeciesConcentrations_j;
      mobileAggregatePrimarySpeciesConcentrations[i] += params.stoichiometricMatrix( j, i+numSecondarySpecies ) * secondarySpeciesConcentrations_j * params.mobileSecondarySpeciesFlag( j );
      for( int k=0; k<numPrimarySpecies; ++k )
      {
        REAL_TYPE const dSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentration = logmath::dWrtLogConst() * params.stoichiometricMatrix( j, k+numSecondarySpecies ) * secondarySpeciesConcentrations_j;
        dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, k ) += params.stoichiometricMatrix( j, i+numSecondarySpecies ) *
                                                                                                      dSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentration;

        dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, k ) += params.stoichiometricMatrix( j, i+numSecondarySpecies ) *
                                                                                                            dSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentration *
                                                                                                            params.mobileSecondarySpeciesFlag( j );
      }
    }
  }
}

} // namespace massActions
} // namespace hpcReact
