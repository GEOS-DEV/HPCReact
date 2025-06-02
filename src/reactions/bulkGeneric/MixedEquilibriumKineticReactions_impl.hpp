#pragma once

#include "common/constants.hpp"
#include "common/CArrayWrapper.hpp"
#include "SpeciesUtilities.hpp"


/** @file MixedEquilibriumKineticReactions_impl.hpp
 *  @brief Header file for the MixedEquilibriumKineticReactions implementation.
 *  @author HPC-REACT Team
 *  @date 2025
 */

namespace hpcReact
{
namespace bulkGeneric
{

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          bool LOGE_CONCENTRATION >  
template< typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_PRIMARY,
          typename ARRAY_1D_SECONDARY,
          typename ARRAY_1D_KINETIC,
          typename ARRAY_2D_PRIMARY,
          typename ARRAY_2D_KINETIC >
HPCREACT_HOST_DEVICE inline void
MixedEquilibriumKineticReactions< REAL_TYPE,
                                  INT_TYPE,
                                  INDEX_TYPE,
                                  LOGE_CONCENTRATION 
                                >::updateMixedSystem_impl( RealType const & temperature,
                                                           PARAMS_DATA const & params,
                                                           ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                           ARRAY_1D_SECONDARY & logSecondarySpeciesConcentrations,
                                                           ARRAY_1D_PRIMARY & aggregatePrimarySpeciesConcentrations,
                                                           ARRAY_2D_PRIMARY & dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                                                           ARRAY_1D_KINETIC & reactionRates,
                                                           ARRAY_2D_KINETIC & dReactionRates_dLogPrimarySpeciesConcentrations,
                                                           ARRAY_1D_PRIMARY & aggregateSpeciesRates,
                                                           ARRAY_2D_PRIMARY & dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations  )
  {
    
    // 1. Compute new aggregate species from primary species
    calculateAggregatePrimaryConcentrationsWrtLogC< REAL_TYPE,
                                                    INT_TYPE,
                                                    INDEX_TYPE >( params.equilibriumReactionsParameters(), 
                                                                  logPrimarySpeciesConcentrations,
                                                                  logSecondarySpeciesConcentrations,
                                                                  aggregatePrimarySpeciesConcentrations,
                                                                  dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations );

    // 2. Compute the reaction rates for all kinetic reactions
    computeReactionRates( temperature,
                          params,
                          logPrimarySpeciesConcentrations,
                          logSecondarySpeciesConcentrations,
                          reactionRates,
                          dReactionRates_dLogPrimarySpeciesConcentrations );
    

    // 3. Compute aggregate species rates
    computeAggregateSpeciesRates( params,
                                  logPrimarySpeciesConcentrations,
                                  reactionRates,
                                  dReactionRates_dLogPrimarySpeciesConcentrations,
                                  aggregateSpeciesRates,
                                  dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations );
    
    
  }
  
  template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          bool LOGE_CONCENTRATION >  
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST2,
            typename ARRAY_1D,
            typename ARRAY_2D >
  HPCREACT_HOST_DEVICE inline void
  MixedEquilibriumKineticReactions< REAL_TYPE,
                                    INT_TYPE,
                                    INDEX_TYPE,
                                    LOGE_CONCENTRATION 
                                     >::computeReactionRates_impl( RealType const & temperature, 
                                                                   PARAMS_DATA const & params,
                                                                   ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                                   ARRAY_1D_TO_CONST2 const & logSecondarySpeciesConcentrations,
                                                                   ARRAY_1D & reactionRates,
                                                                   ARRAY_2D & dReactionRates_dLogPrimarySpeciesConcentrations )

  {
    constexpr IntType numSpecies          = PARAMS_DATA::numSpecies();
    constexpr IntType numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
    constexpr IntType numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();
    constexpr IntType numKineticReactions = PARAMS_DATA::numKineticReactions();

    RealType logSpeciesConcentration[numSpecies] {};
    for ( IntType i = 0; i < numSecondarySpecies; ++i )
    {
      logSpeciesConcentration[i] = logSecondarySpeciesConcentrations[i];
    }
    for ( IntType i = numSecondarySpecies; i < numSpecies; ++i )
    {
      logSpeciesConcentration[i] = logPrimarySpeciesConcentrations[i];
    }
    
    CArrayWrapper< RealType, numKineticReactions, numSpecies > reactionRatesDerivatives;

    kineticReactions::computeReactionRates( temperature,
                                            params.kineticReactionsParameters(),
                                            logSpeciesConcentration,
                                            reactionRates,
                                            reactionRatesDerivatives );

    // Compute the reaction rates derivatives w.r.t. log primary species concentrations
    for( IntType i = 0; i < numKineticReactions; ++i )
    {
      for( IntType j = 0; j < numPrimarySpecies; ++j )
      { 
        dReactionRates_dLogPrimarySpeciesConcentrations( i, j ) = reactionRatesDerivatives( i, j + numSecondarySpecies );

        for( IntType k = 0; k < numSecondarySpecies; ++k )
        {
          RealType const dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations = params.stoichiometricMatrix( k, j + numSecondarySpecies );
          
          dReactionRates_dLogPrimarySpeciesConcentrations( i, j ) += 
          reactionRatesDerivatives( i, k ) * dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations;
         }
      }
    }
  }                                        

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          bool LOGE_CONCENTRATION >  
template< typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D,
          bool CALCULATE_DERIVATIVES >
HPCREACT_HOST_DEVICE inline void
MixedEquilibriumKineticReactions< REAL_TYPE,
                                  INT_TYPE,
                                  INDEX_TYPE,
                                  LOGE_CONCENTRATION 
                                >::computeAggregateSpeciesRates_impl( PARAMS_DATA const & params,
                                                                      ARRAY_1D_TO_CONST const & speciesConcentration,
                                                                      ARRAY_1D_TO_CONST2 const & reactionRates,
                                                                      ARRAY_2D_TO_CONST const & reactionRatesDerivatives, 
                                                                      ARRAY_1D & aggregatesRates,
                                                                      ARRAY_2D & aggregatesRatesDerivatives )
  {
    HPCREACT_UNUSED_VAR( speciesConcentration );
    // constexpr IntType numSpecies = PARAMS_DATA::numSpecies();
    constexpr IntType numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
    constexpr IntType numKineticReactions = PARAMS_DATA::numKineticReactions();
    constexpr IntType numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

    for( IntType i = 0; i < numPrimarySpecies; ++i )
    {
      aggregatesRates[i] = 0.0;
      if constexpr( CALCULATE_DERIVATIVES )
      {
        for( IntType j = 0; j < numPrimarySpecies; ++j )
        {
          aggregatesRatesDerivatives( i, j ) = 0.0;
       }
      }
      for( IntType r=0; r<numKineticReactions; ++r )
      {
        RealType const s_ir = params.kineticReactionsParameters().stoichiometricMatrix( r, i+numSecondarySpecies );
        aggregatesRates[i] += s_ir * reactionRates[r];
        if constexpr( CALCULATE_DERIVATIVES )
        {
          for( IntType j = 0; j < numPrimarySpecies; ++j )
          {
            aggregatesRatesDerivatives( i, j ) += s_ir * reactionRatesDerivatives( r, j );
          }
        }
      } 
    }

  }

} // namespace bulkGeneric

} // namespace hpcReact

#include "common/macrosCleanup.hpp"
