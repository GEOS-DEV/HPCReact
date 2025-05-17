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
                     ARRAY_2D & dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                     ARRAY_1D & reactionRates,
                     ARRAY_2D & dReactionRates_dLogPrimarySpeciesConcentrations,
                     ARRAY_1D & aggregateSpeciesRates,
                     ARRAY_2D & dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations )
  {
    updateMixedSystem_impl( temperature,
                            params,
                            logPrimarySpeciesConcentrations,
                            logSecondarySpeciesConcentrations,
                            aggregatePrimarySpeciesConcentrations,
                            dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
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
                        ARRAY_1D & reactionRates,
                        ARRAY_2D & dReactionRates_dPrimarySpeciesConcentrations )

  {
    computeReactionRates_impl( temperature,
                               params,
                               logPrimarySpeciesConcentrations,
                               logSecondarySpeciesConcentrations,
                               reactionRates,
                               dReactionRates_dPrimarySpeciesConcentrations );
  }                                        

  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_2D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeAggregateSpeciesRates( PARAMS_DATA const & params,
                                ARRAY_1D_TO_CONST const & speciesConcentration,
                                ARRAY_1D_TO_CONST const & reactionRates,
                                ARRAY_2D_TO_CONST const & reactionRatesDerivatives, 
                                ARRAY_1D & aggregatesRates,
                                ARRAY_2D & aggregatesRatesDerivatives )
  {
    computeAggregateSpeciesRates_impl< true >( params,
                                               speciesConcentration,
                                               reactionRates,
                                               reactionRatesDerivatives,
                                               aggregatesRates,
                                               aggregatesRatesDerivatives );
  }

  private:

  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  updateMixedSystem_impl( RealType const & temperature,
                          PARAMS_DATA const & params,
                          ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                          ARRAY_1D & logSecondarySpeciesConcentrations,
                          ARRAY_1D & aggregatePrimarySpeciesConcentrations,
                          ARRAY_2D & dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                          ARRAY_1D & reactionRates,
                          ARRAY_2D & dReactionRates_dLogPrimarySpeciesConcentrations,
                          ARRAY_1D & aggregateSpeciesRates,
                          ARRAY_2D & dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations );

  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeReactionRates_impl( RealType const & temperature,
                        PARAMS_DATA const & params,
                        ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                        ARRAY_1D_TO_CONST const & logSecondarySpeciesConcentrations,
                        ARRAY_1D & reactionRates,
                        ARRAY_2D & dReactionRates_dPrimarySpeciesConcentrations );



  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_2D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D,
            bool CALCULATE_DERIVATIVES >
  static HPCREACT_HOST_DEVICE void
  computeAggregateSpeciesRates_impl( PARAMS_DATA const & params,
                                     ARRAY_1D_TO_CONST const & speciesConcentration,
                                     ARRAY_1D_TO_CONST const & reactionRates,
                                     ARRAY_2D_TO_CONST const & reactionRatesDerivatives, 
                                     ARRAY_1D & aggregatesRates,
                                     ARRAY_2D & aggregatesRatesDerivatives );


};

} // namespace bulkGeneric
} // namespace hpcReact

#include "MixedEquilibriumKineticReactions_impl.hpp"
#include "common/macrosCleanup.hpp"
