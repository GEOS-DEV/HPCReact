#pragma once

#include "common/macros.hpp"
#include "KineticReactions.hpp"

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
  
  /// Type alias for the Kinetic reactions type used in the class.
  using kineticReactions = KineticReactions< REAL_TYPE, INT_TYPE, INDEX_TYPE, LOGE_CONCENTRATION >;
  
  /**
   * @brief Update a mixed chemical system by computing secondary species concentrations,
   * aggregate primary species concentrations, and reaction rates.
   *
   * @tparam PARAMS_DATA Struct providing all parameter access (stoichiometry, rate constants, etc.)
   * @tparam ARRAY_1D_TO_CONST Read-only 1D array type for primary log-concentrations
   * @tparam ARRAY_1D_PRIMARY Mutable 1D array type for primary species outputs
   * @tparam ARRAY_1D_SECONDARY Mutable 1D array type for secondary log-concentrations
   * @tparam ARRAY_1D_KINETIC Mutable 1D array type for reaction rates
   * @tparam ARRAY_2D_PRIMARY Mutable 2D array type for primary derivatives
   * @tparam ARRAY_2D_KINETIC Mutable 2D array type for reaction rate derivatives
   *
   * @param temperature Temperature of the system (in Kelvin)
   * @param params Parameter object for stoichiometry, rates, etc.
   * @param logPrimarySpeciesConcentrations Log of primary species concentrations
   * @param surfaceArea surface Aread for kinetic reactions
   * @param logSecondarySpeciesConcentrations Output log concentrations for secondary species
   * @param aggregatePrimarySpeciesConcentrations Output aggregate concentrations (per primary)
   * @param dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations Derivatives of aggregate concentrations w.r.t. log primary
   * @param reactionRates Output vector of kinetic reaction rates
   * @param dReactionRates_dLogPrimarySpeciesConcentrations Derivatives of reaction rates w.r.t. log primary species
   * @param aggregateSpeciesRates Output net source/sink for each primary species
   * @param dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations Derivatives of aggregate source terms
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST_KINETIC,
            typename ARRAY_1D_PRIMARY,
            typename ARRAY_1D_SECONDARY,
            typename ARRAY_1D_KINETIC,
            typename ARRAY_2D_PRIMARY,
            typename ARRAY_2D_KINETIC >
  static HPCREACT_HOST_DEVICE inline void
  updateMixedSystem( RealType const & temperature,
                     PARAMS_DATA const & params,
                     ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                     ARRAY_1D_TO_CONST_KINETIC const & surfaceArea,
                     ARRAY_1D_SECONDARY & logSecondarySpeciesConcentrations,
                     ARRAY_1D_PRIMARY & aggregatePrimarySpeciesConcentrations,
                     ARRAY_2D_PRIMARY & dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                     ARRAY_1D_KINETIC & reactionRates,
                     ARRAY_2D_KINETIC & dReactionRates_dLogPrimarySpeciesConcentrations,
                     ARRAY_1D_PRIMARY & aggregateSpeciesRates,
                     ARRAY_2D_PRIMARY & dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations )
  {
    updateMixedSystem_impl( temperature,
                            params,
                            logPrimarySpeciesConcentrations,
                            surfaceArea,
                            logSecondarySpeciesConcentrations,
                            aggregatePrimarySpeciesConcentrations,
                            dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                            reactionRates,
                            dReactionRates_dLogPrimarySpeciesConcentrations,
                            aggregateSpeciesRates,
                            dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations );
  }
  
 /**
  * @brief Compute reaction rates and their derivatives.
  *
  * @tparam PARAMS_DATA Struct providing reaction parameters
  * @tparam ARRAY_1D_TO_CONST Read-only array of primary species (log-space)
  * @tparam ARRAY_1D_TO_CONST2 Read-only array of secondary species (log-space)
  * @tparam ARRAY_1D Output array type for reaction rates
  * @tparam ARRAY_2D Output array type for reaction rate derivatives
  *
  * @param temperature Temperature in Kelvin
  * @param params Parameter data for the reaction system
  * @param logPrimarySpeciesConcentrations Log concentrations of primary species
  * @param logSecondarySpeciesConcentrations Log concentrations of secondary species
  * @param surfaceArea Surface area for kinetic reactions
  * @param reactionRates Output reaction rates for each kinetic reaction
  * @param dReactionRates_dLogPrimarySpeciesConcentrations Derivatives of reaction rates w.r.t. log primary species
  */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST2,
            typename ARRAY_1D_TO_CONST_KINETIC,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                        ARRAY_1D_TO_CONST2 const & logSecondarySpeciesConcentrations,
                        ARRAY_1D_TO_CONST_KINETIC const & surfaceArea,
                        ARRAY_1D & reactionRates,
                        ARRAY_2D & dReactionRates_dLogPrimarySpeciesConcentrations )

  {
    computeReactionRates_impl( temperature,
                               params,
                               logPrimarySpeciesConcentrations,
                               logSecondarySpeciesConcentrations,
                               surfaceArea,
                               reactionRates,
                               dReactionRates_dLogPrimarySpeciesConcentrations );
  }                                        
  
  /**
   * @brief Compute net reaction rate for each primary species by aggregating contributions from all reactions.
   *
   * @tparam PARAMS_DATA Struct with stoichiometry and mappings
   * @tparam ARRAY_1D_TO_CONST Array type for primary species concentrations
   * @tparam ARRAY_1D_TO_CONST2 Array type for reaction rates
   * @tparam ARRAY_2D_TO_CONST Array type for reaction rate derivatives
   * @tparam ARRAY_1D Output type for net species rates
   * @tparam ARRAY_2D Output type for rate derivatives
   *
   * @param params Reaction parameters
   * @param speciesConcentration Current concentrations of primary species
   * @param reactionRates Computed reaction rates
   * @param reactionRatesDerivatives Derivatives of reaction rates w.r.t. log concentrations
   * @param aggregatesRates Output: net rate for each primary species
   * @param aggregatesRatesDerivatives Output: derivative of net rates w.r.t. log concentrations
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST2,
            typename ARRAY_2D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeAggregateSpeciesRates( PARAMS_DATA const & params,
                                ARRAY_1D_TO_CONST const & speciesConcentration,
                                ARRAY_1D_TO_CONST2 const & reactionRates,
                                ARRAY_2D_TO_CONST const & reactionRatesDerivatives, 
                                ARRAY_1D & aggregatesRates,
                                ARRAY_2D & aggregatesRatesDerivatives )
  {
    computeAggregateSpeciesRates_impl<  PARAMS_DATA, 
                                        ARRAY_1D_TO_CONST,
                                        ARRAY_1D_TO_CONST2,
                                        ARRAY_2D_TO_CONST,
                                        ARRAY_1D,
                                        ARRAY_2D,
                                        true >( params,
                                               speciesConcentration,
                                               reactionRates,
                                               reactionRatesDerivatives,
                                               aggregatesRates,
                                               aggregatesRatesDerivatives );
  }

  private:
  /**
   * @brief Internal implementation of updateMixedSystem with template-dispatched logic.
   *
   * @details Called by the public `updateMixedSystem` function. Handles the complete chain:
   *          secondary speciation, aggregation, reaction rate evaluation, and net source terms.
   * @tparam PARAMS_DATA Struct providing all parameter access (stoichiometry, rate constants, etc.)
   * @tparam ARRAY_1D_TO_CONST Read-only 1D array type for primary log-concentrations
   * @tparam ARRAY_1D_PRIMARY Mutable 1D array type for primary species outputs
   * @tparam ARRAY_1D_SECONDARY Mutable 1D array type for secondary log-concentrations
   * @tparam ARRAY_1D_KINETIC Mutable 1D array type for reaction rates
   * @tparam ARRAY_2D_PRIMARY Mutable 2D array type for primary derivatives
   * @tparam ARRAY_2D_KINETIC Mutable 2D array type for reaction rate derivatives
   *
   * @param temperature Temperature of the system (in Kelvin)
   * @param params Parameter object for stoichiometry, rates, etc.
   * @param logPrimarySpeciesConcentrations Log of primary species concentrations
   * @param logSecondarySpeciesConcentrations Output log concentrations for secondary species
   * @param aggregatePrimarySpeciesConcentrations Output aggregate concentrations (per primary)
   * @param dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations Derivatives of aggregate concentrations w.r.t. log primary
   * @param reactionRates Output vector of kinetic reaction rates
   * @param dReactionRates_dLogPrimarySpeciesConcentrations Derivatives of reaction rates w.r.t. log primary species
   * @param aggregateSpeciesRates Output net source/sink for each primary species
   * @param dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations Derivatives of aggregate source terms
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST_KINETIC,
            typename ARRAY_1D_PRIMARY,
            typename ARRAY_1D_SECONDARY,
            typename ARRAY_1D_KINETIC, 
            typename ARRAY_2D_PRIMARY,
            typename ARRAY_2D_KINETIC >
  static HPCREACT_HOST_DEVICE void
  updateMixedSystem_impl( RealType const & temperature,
                          PARAMS_DATA const & params,
                          ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                          ARRAY_1D_TO_CONST_KINETIC const & surfaceArea,
                          ARRAY_1D_SECONDARY & logSecondarySpeciesConcentrations,
                          ARRAY_1D_PRIMARY & aggregatePrimarySpeciesConcentrations,
                          ARRAY_2D_PRIMARY & dAggregatePrimarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                          ARRAY_1D_KINETIC & reactionRates,
                          ARRAY_2D_KINETIC & dReactionRates_dLogPrimarySpeciesConcentrations,
                          ARRAY_1D_PRIMARY & aggregateSpeciesRates,
                          ARRAY_2D_PRIMARY & dAggregateSpeciesRates_dLogPrimarySpeciesConcentrations );
  /**
   * @brief Internal implementation of computeReactionRates.
   *
   * @details Handles kinetic rate law evaluation for forward and reverse reactions.
   * @param temperature Temperature in Kelvin
   * @param params Parameter data for the reaction system
   * @param logPrimarySpeciesConcentrations Log concentrations of primary species
   * @param logSecondarySpeciesConcentrations Log concentrations of secondary species
   * @param reactionRates Output reaction rates for each kinetic reaction
   * @param dReactionRates_dLogPrimarySpeciesConcentrations Derivatives of reaction rates w.r.t. log primary species
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST2,
            typename ARRAY_1D_TO_CONST_KINETIC,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeReactionRates_impl( RealType const & temperature,
                        PARAMS_DATA const & params,
                        ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                        ARRAY_1D_TO_CONST2 const & logSecondarySpeciesConcentrations,
                        ARRAY_1D_TO_CONST_KINETIC const & surfaceArea,
                        ARRAY_1D & reactionRates,
                        ARRAY_2D & dReactionRates_dLogPrimarySpeciesConcentrations );


  /**
   * @brief Internal implementation of computeAggregateSpeciesRates.
   *
   * @tparam CALCULATE_DERIVATIVES Whether to compute Jacobian derivatives
   * @tparam ARRAY_1D_TO_CONST Array type for primary species concentrations
   * @tparam ARRAY_1D_TO_CONST2 Array type for reaction rates
   * @tparam ARRAY_2D_TO_CONST Array type for reaction rate derivatives
   * @tparam ARRAY_1D Output type for net species rates
   * @tparam ARRAY_2D Output type for rate derivatives
   *
   * @param params Reaction parameters
   * @param speciesConcentration Current concentrations of primary species
   * @param reactionRates Computed reaction rates
   * @param reactionRatesDerivatives Derivatives of reaction rates w.r.t. log concentrations
   * @param aggregatesRates Output: net rate for each primary species
   * @param aggregatesRatesDerivatives Output: derivative of net rates w.r.t. log concentrations
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST2,
            typename ARRAY_2D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D,
            bool CALCULATE_DERIVATIVES >
  static HPCREACT_HOST_DEVICE void
  computeAggregateSpeciesRates_impl( PARAMS_DATA const & params,
                                     ARRAY_1D_TO_CONST const & speciesConcentration,
                                     ARRAY_1D_TO_CONST2 const & reactionRates,
                                     ARRAY_2D_TO_CONST const & reactionRatesDerivatives, 
                                     ARRAY_1D & aggregatesRates,
                                     ARRAY_2D & aggregatesRatesDerivatives );


};

} // namespace bulkGeneric
} // namespace hpcReact

#include "MixedEquilibriumKineticReactions_impl.hpp"
#include "common/macrosCleanup.hpp"
