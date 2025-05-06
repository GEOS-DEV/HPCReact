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
          bool LOGE_CONCENTRATION >
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
                     ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations,
                     ARRAY_1D & reactionRates,
                     ARRAY_2D & reactionRatesDerivatives )
  {
    constexpr int numSpecies = PARAMS_DATA::numSpecies;
    constexpr int numSecondarySpecies = PARAMS_DATA::numReactions;
    constexpr int numPrimarySpecies = numSpecies - numSecondarySpecies;

    // 1. Compute new aggregate species from primary species
    calculateAggregatePrimaryConcentrationsWrtLogC< REAL_TYPE,
                                                    INT_TYPE,
                                                    INDEX_TYPE >( params, 
                                                                  logPrimarySpeciesConcentrations,
                                                                  logSecondarySpeciesConcentrations,
                                                                  aggregatePrimarySpeciesConcentrations,
                                                                  dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )

    // 2. Compute the reaction rates for all kinetic reactions
    kineticReactions::computeReactionRates( temperature,
                                            params,
                                            speciesConcentration,
                                            reactionRates,
                                            reactionRatesDerivatives );

    // 3. Compute aggregate species rates
    
    
  }

private:

};

} // namespace bulkGeneric
} // namespace hpcReact

#include "MixedEquilibriumKineticReactions_impl.hpp"
#include "common/macrosCleanup.hpp"
