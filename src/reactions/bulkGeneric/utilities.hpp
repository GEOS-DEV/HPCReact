#pragma once

#include "common/macros.hpp"
#include <math.h>

namespace hpcReact
{
namespace bulkGeneric
{

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
void calculateAggregatePrimaryConcentrationsWrtLogC( PARAMS_DATA const & params,
                                                     ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                     ARRAY_1D & aggregatePrimarySpeciesConcentrations,
                                                     ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  constexpr int numSpecies = PARAMS_DATA::numSpecies;
  constexpr int numSecondarySpecies = PARAMS_DATA::numReactions;
  constexpr int numPrimarySpecies = numSpecies - numSecondarySpecies;

  REAL_TYPE logSecondarySpeciesConcentrations[numSecondarySpecies] = {0};

  for( int j=0; j<numSecondarySpecies; ++j )
  {
    REAL_TYPE const gamma = 1;
    logSecondarySpeciesConcentrations[j] = -log( params.equilibriumConstant( j ) ) - log( gamma );
    for( int k=0; k<numPrimarySpecies; ++k )
    {
      logSecondarySpeciesConcentrations[j] += params.stoichiometricMatrix( j, k+numSecondarySpecies ) * ( logPrimarySpeciesConcentrations[k] + log( gamma ) );
    }
  }

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
void calculateAggregateBasedResidualAndJacobianWrtLogC( )

} // namespace bulkGeneric
} // namespace hpcReact
