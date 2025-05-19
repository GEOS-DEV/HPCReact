#pragma once

#include "common/macros.hpp"
#include <math.h>
#include <functional>
#include <iostream>

namespace hpcReact
{
namespace bulkGeneric
{

namespace utilities_impl
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
  constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies(); 

  for( int j=0; j<numSecondarySpecies; ++j )
  {
    REAL_TYPE const gamma = 1;
    logSecondarySpeciesConcentrations[j] = -log( params.equilibriumConstant( j ) ) - log( gamma );
    for( int k=0; k<numPrimarySpecies; ++k )
    {
      logSecondarySpeciesConcentrations[j] += params.stoichiometricMatrix( j, k+numSecondarySpecies ) * ( logPrimarySpeciesConcentrations[k] + log( gamma ) );
      derivativeFunc( j, k, params.stoichiometricMatrix( j, k+numSecondarySpecies ) );
    }
  }
}

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
  utilities_impl::calculateLogSecondarySpeciesConcentration< REAL_TYPE,
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
  utilities_impl::calculateLogSecondarySpeciesConcentration< REAL_TYPE, INT_TYPE, INDEX_TYPE >( params,
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
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateAggregatePrimaryConcentrationsWrtLogC( PARAMS_DATA & params,
                                                     ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                     ARRAY_1D & logSecondarySpeciesConcentrations,
                                                     ARRAY_1D & aggregatePrimarySpeciesConcentrations,
                                                     ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();
  constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();


  calculateLogSecondarySpeciesConcentration< REAL_TYPE,
                                             INT_TYPE,
                                             INDEX_TYPE >( params,
                                                           logPrimarySpeciesConcentrations,
                                                           logSecondarySpeciesConcentrations );

  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    REAL_TYPE const speciesConcentration_i = exp( logPrimarySpeciesConcentrations[i] );
    aggregatePrimarySpeciesConcentrations[i] = speciesConcentration_i;
    dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, i ) = speciesConcentration_i;
    for( int j = 0; j < numSecondarySpecies; ++j )
    {
      REAL_TYPE const secondarySpeciesConcentrations_j = exp( logSecondarySpeciesConcentrations[j] );
      aggregatePrimarySpeciesConcentrations[i] += params.stoichiometricMatrix( j, i+numSecondarySpecies ) * secondarySpeciesConcentrations_j;
      for( int k=0; k<numPrimarySpecies; ++k )
      {
        REAL_TYPE const dSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentration = params.stoichiometricMatrix( j, k+numSecondarySpecies ) * secondarySpeciesConcentrations_j;
        dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, k ) += params.stoichiometricMatrix( j,
                                                                                                                                   i+numSecondarySpecies ) *
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
  constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  REAL_TYPE logSecondarySpeciesConcentrations[numSecondarySpecies] = {0};

  calculateLogSecondarySpeciesConcentration< REAL_TYPE,
                                             INT_TYPE,
                                             INDEX_TYPE >( params,
                                                           logPrimarySpeciesConcentrations,
                                                           logSecondarySpeciesConcentrations );

  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    REAL_TYPE const speciesConcentration_i = exp( logPrimarySpeciesConcentrations[i] );
    aggregatePrimarySpeciesConcentrations[i] = speciesConcentration_i;
    dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, i ) = speciesConcentration_i;
    for( int j = 0; j < numSecondarySpecies; ++j )
    {
      REAL_TYPE const secondarySpeciesConcentrations_j = exp( logSecondarySpeciesConcentrations[j] );
      aggregatePrimarySpeciesConcentrations[i] += params.stoichiometricMatrix( j, i+numSecondarySpecies ) * secondarySpeciesConcentrations_j;
      for( int k=0; k<numPrimarySpecies; ++k )
      {
        REAL_TYPE const dSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentration = params.stoichiometricMatrix( j, k+numSecondarySpecies ) * secondarySpeciesConcentrations_j;
        dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, k ) += params.stoichiometricMatrix( j,
                                                                                                                                   i+numSecondarySpecies ) *
                                                                                                      dSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentration;
      }
    }
  }
}

} // namespace bulkGeneric
} // namespace hpcReact
