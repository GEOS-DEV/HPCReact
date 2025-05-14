#pragma once

#include "common/macros.hpp"

/** @file MixedEquilibriumKineticReactions.hpp
 *  @brief Header file for the MixedEquilibriumKineticReactions class.
 *  @author HPC-REACT Team
 *  @date 2025
 */

namespace hpcReact
{
namespace bulkGeneric
{

/**
 * @brief Class for computing reaction rates and species rates for a given set of reactions.
 * @tparam REAL_TYPE The type of the real numbers used in the class.
 * @tparam INT_TYPE The type of the integer numbers used in the class.
 * @tparam INDEX_TYPE The type of the index used in the class.
 * @tparam LOGE_CONCENTRATION Whether to use logarithm of concentration for the calculations.
 * @details
 *   This class provides the ablity to compute kinetic reactions.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          bool LOGE_CONCENTRATION = true >
class MixedEquilibriumKineticReactions
{
public:

  /// Type alias for the real type used in the class.
  using RealType = REAL_TYPE;

  /// Type alias for the integer type used in the class.
  using IntType = INT_TYPE;

  /// Type alias for the index type used in the class.
  using IndexType = INDEX_TYPE;

  using kineticReactions = KineticReactions< REAL_TYPE, INT_TYPE, INDEX_TYPE, LOGE_CONCENTRATION >;

  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  updateMixedSystem( RealType const & temperature,
                     PARAMS_DATA const & params,
                     ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                     ARRAY_1D & logSecondarySpeciesConcentrations,
                     ARRAY_1D & aggregatePrimarySpeciesConcentrations,
                     ARRAY_2D & dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                     ARRAY_1D & reactionRates,
                     ARRAY_2D & dReactionRates_dLogPrimarySpeciesConcentrations,
                     ARRAY_1D & aggregateSpeciesRates,
                     ARRAY_2D & dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations )
  {
    constexpr IntType numSpecies = PARAMS_DATA::numSpecies;
    constexpr IntType numSecondarySpecies = PARAMS_DATA::numReactions;
    constexpr IntType numPrimarySpecies = numSpecies - numSecondarySpecies;

    // 1. Compute new aggregate species from primary species
    calculateAggregatePrimaryConcentrationsWrtLogC< REAL_TYPE,
                                                    INT_TYPE,
                                                    INDEX_TYPE >( params, 
                                                                  logPrimarySpeciesConcentrations,
                                                                  logSecondarySpeciesConcentrations,
                                                                  aggregatePrimarySpeciesConcentrations,
                                                                  dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )

    // 2. Compute the reaction rates for all kinetic reactions
    computeReactionRates( temperature,
                          params,
                          logPrimarySpeciesConcentrations,
                          logSecondarySpeciesConcentrations,
                          aggregatePrimarySpeciesConcentrations,
                          dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                          reactionRates,
                          dReactionRates_dLogPrimarySpeciesConcentrations );
    

    // 3. Compute aggregate species rates
    computeAggregateSpeciesRates( params,
                                  logPrimarySpeciesConcentrations
                                  reactionRates,
                                  dReactionRates_dLogPrimarySpeciesConcentrations,
                                  aggregateSpeciesRates,
                                  dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations );
    
    
  }
  
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                        ARRAY_1D_TO_CONST const & logSecondarySpeciesConcentrations,
                        ARRAY_2D_TO_CONST const & dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                        ARRAY_1D & reactionRates,
                        ARRAY_2D & dReactionRates_dPrimarySpeciesConcentrations )

  {
    constexpr IntType numSpecies = PARAMS_DATA::numSpecies;
    constexpr IntType numSecondarySpecies = PARAMS_DATA::numReactions;
    constexpr IntType numPrimarySpecies = numSpecies - numSecondarySpecies;

    RealType logSpeciesConcentration[numSpecies] {};
    for ( IntType i = 0; i < numPrimarySpecies; ++i )
    {
      speciesConcentration[i] = logPrimarySpeciesConcentrations[i];
    }
    for ( IntType i = 0; i < numSecondarySpecies; ++i )
    {
      logSpeciesConcentration[i+numPrimarySpecies] = logSecondarySpeciesConcentrations[i];
    }
    
    CArrayWrapper< RealType, PARAMS_DATA::numReactions, numSpecies > reactionRatesDerivatives;
    kineticReactions::computeReactionRates( temperature,
                                            params,
                                            logSpeciesConcentration,
                                            reactionRates,
                                            reactionRatesDerivatives );

    // Compute the reaction rates derivatives w.r.t. log primary species concentrations
    for( IntType i = 0; i < numKineticReactions; ++i )
    {
      for( IntType j = 0; j < numPrimarySpecies; ++j )
      {
        dReactionRates_dLogPrimarySpeciesConcentrations( i, j ) = reactionRatesDerivatives( i, j );
        for( IntType k = 0; k < numSecondarySpecies; ++k )
        {
          RealType const dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentration = params.stoichiometricMatrix( k, j+numSecondarySpecies );
          
          dReactionRates_dLogPrimarySpeciesConcentrations( i, j ) += 
          reactionRatesDerivatives( i, numPrimarySpecies + k ) * dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations( k, j );
        }
      }
    }
  }                                        

  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_2D_TO_CONST
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeAggregateSpeciesRates( PARAMS_DATA const & params,
                                ARRAY_1D_TO_CONST const & speciesConcentration,
                                ARRAY_1D_TO_CONST const & reactionRates,
                                ARRAY_2D_TO_CONST const & reactionRatesDerivatives 
                                ARRAY_1D & aggregatesRates,
                                ARRAY_2D & aggregatesRatesDerivatives )
  {
    constexpr IntType numSpecies = PARAMS_DATA::numSpecies;
    constexpr IntType numSecondarySpecies = PARAMS_DATA::numReactions;
    constexpr IntType numPrimarySpecies = numSpecies - numSecondarySpecies;

    for( IntType i = 0; i < numPrimarySpecies; ++i )
    {
      speciesRates[i] = 0.0;
      if constexpr( CALCULATE_DERIVATIVES )
      {
        for( IntType j = 0; j < numPrimarySpecies; ++j )
        {
          speciesRatesDerivatives( i, j ) = 0.0;
       }
      }
      for( IntType r=0; r<PARAMS_DATA::numReactions; ++r )
      {
        RealType const s_ir = params.stoichiometricMatrix( r, i );
        speciesRates[i] += s_ir * reactionRates[r];
        if constexpr( CALCULATE_DERIVATIVES )
        {
          for( IntType j = 0; j < numPrimarySpecies; ++j )
          {
            speciesRatesDerivatives( i, j ) += s_ir * reactionRatesDerivatives( r, j );
          }
        }
      } 
    }

  }

private:

};

} // namespace bulkGeneric
} // namespace hpcReact

#include "MixedEquilibriumKineticReactions_impl.hpp"
#include "common/macrosCleanup.hpp"
