#pragma once

#include "common/macros.hpp"

/** @file KineticReactions.hpp
 *  @brief Header file for the KineticReactions class.
 *  @author HPC-REACT Team
 *  @date 2023
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
class KineticReactions
{
public:

  /// Type alias for the real type used in the class.
  using RealType = REAL_TYPE;

  /// Type alias for the integer type used in the class.
  using IntType = INT_TYPE;

  /// Type alias for the index type used in the class.
  using IndexType = INDEX_TYPE;

  /**
   * @copydoc KineticReactions::computeReactionRates_impl()
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        ARRAY_1D_TO_CONST const & speciesConcentration,
                        ARRAY_1D & reactionRates,
                        ARRAY_2D & reactionRatesDerivatives )
  {
    computeReactionRates_impl< PARAMS_DATA, true >( temperature,
                                                    params,
                                                    speciesConcentration,
                                                    reactionRates,
                                                    reactionRatesDerivatives );
  }

  /**
   * @brief Compute the reaction rates for a given set of species concentrations.
   * @tparam PARAMS_DATA The type of the parameters data.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @tparam ARRAY_1D The type of the array of reaction rates.
   * @param temperature The temperature of the system.
   * @param params The parameters data.
   * @param speciesConcentration The array of species concentrations.
   * @param reactionRates The array of reaction rates.
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        ARRAY_1D_TO_CONST const & speciesConcentration,
                        ARRAY_1D & reactionRates )
  {
    REAL_TYPE reactionRatesDerivatives[PARAMS_DATA::numReactions][PARAMS_DATA::numSpecies] = { {0.0} };
    computeReactionRates_impl< PARAMS_DATA, false >( temperature,
                                                     params,
                                                     speciesConcentration,
                                                     reactionRates,
                                                     reactionRatesDerivatives );
  }


  /**
   * @copydoc KineticReactions::computeSpeciesRates_impl()
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeSpeciesRates( RealType const & temperature,
                       PARAMS_DATA const & params,
                       ARRAY_1D_TO_CONST const & speciesConcentration,
                       ARRAY_1D & speciesRates,
                       ARRAY_2D & speciesRatesDerivatives )
  {
    computeSpeciesRates_impl< PARAMS_DATA, true >( temperature,
                                                   params,
                                                   speciesConcentration,
                                                   speciesRates,
                                                   speciesRatesDerivatives );
  }

  /**
   * @brief Compute the species rates for a given set of species concentrations.
   * @tparam PARAMS_DATA The type of the parameters data.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @tparam ARRAY_1D The type of the array of species rates.
   * @param temperature The temperature of the system.
   * @param params The parameters data.
   * @param speciesConcentration The array of species concentrations.
   * @param speciesRates The array of species rates.
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D >
  static HPCREACT_HOST_DEVICE inline void
  computeSpeciesRates( RealType const & temperature,
                       PARAMS_DATA const & params,
                       ARRAY_1D_TO_CONST const & speciesConcentration,
                       ARRAY_1D & speciesRates )
  {
    char speciesRatesDerivatives;
    computeSpeciesRates_impl< PARAMS_DATA, false >( temperature,
                                                    params,
                                                    speciesConcentration,
                                                    speciesRates,
                                                    speciesRatesDerivatives );
  }

  /**
   * @brief execute the time step for a given set of kinetic reactions.
   * @tparam PARAMS_DATA The type of the parameters data.
   * @tparam ARRAY_1D The type of the array of species concentrations.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @tparam ARRAY_2D The type of the array of species rates derivatives.
   * @param dt The time step to be used for the simulation.
   * @param temperature The temperature of the reaction.
   * @param params The parameters data.
   * @param speciesConcentration_n The array of species concentrations at the beginning of the time step.
   * @param speciesConcentration The array of species concentrations at the end of the time step.
   * @param speciesRates The array of species rates.
   * @param speciesRatesDerivatives The array of species rates derivatives.
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  timeStep( RealType const dt,
            RealType const & temperature,
            PARAMS_DATA const & params,
            ARRAY_1D_TO_CONST const & speciesConcentration_n,
            ARRAY_1D & speciesConcentration,
            ARRAY_1D & speciesRates,
            ARRAY_2D & speciesRatesDerivatives );


private:

  /**
   * @brief Compute the reaction rates for a given set of species concentrations.
   * @tparam PARAMS_DATA The type of the parameters data.
   * @tparam CALCULATE_DERIVATIVES Whether to calculate the derivatives of the reaction rates with respect to the species concentrations.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @tparam ARRAY_1D The type of the array of reaction rates.
   * @tparam ARRAY_2D The type of the array of reaction rates derivatives.
   * @param temperature The temperature of the system.
   * @param params The parameters data.
   * @param speciesConcentration The array of species concentrations.
   * @param reactionRates The array of reaction rates.
   * @param reactionRatesDerivatives The array of reaction rates derivatives.
   * @details
   *   This function computes the reaction rates for a given set of reactions.
   *   If CALCULATE_DERIVATIVES is true, it also computes the derivatives of the reaction rates with respect to the species concentrations.
   *
   * The expression for the reaction rate ( \f$ \dot{R}_r \f$ ) for reaction \f$ r \f$  is given by:
   *  \f[
          \dot{R}_r = k^f_{r} \prod_{ \mathclap{ \substack{ i = 1,...,N_s \\ \nu_{ri} < 0 } } } [C_i]^{-\nu_{ri}}
                    - k^r_{r} \prod_{ \mathclap{ \substack{ i = 1,...,N_s \\ \nu_{ri} > 0 } } } [C_i]^{\nu_{ri}}
   *  \f]
   *
   * The expression for the derivative of the reaction rate ( \f$ \frac{d\dot{R}_r}{dC_i} \f$ ) for reaction \f$ r \f$  is given by:
   *  \f[
   *    \frac{\partial \dot{R}_r}{ \partial [C]_i } =
   *      k^f_{r} \left( -\nu_{ri} [C_i]^{-\nu_{ri}-1 } \prod_{ \mathclap{ \substack{ j = 1,...,N_s \\ \nu_{rj} < 0, j \ne i } } }
   * [C_j]^{-\nu_{rj}} \right )
   *    - k^r_{r} \left(  \nu_{ri} [C_i]^{ \nu_{ri}-1 }  \prod_{ \mathclap{ \substack{ j = 1,...,N_s \\ \nu_{rj} > 0, j \ne i } } }
   * [C_j]^{\nu_{rj}} \right )
   *  \f]
   */
  template< typename PARAMS_DATA,
            bool CALCULATE_DERIVATIVES,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeReactionRates_impl( RealType const & temperature,
                             PARAMS_DATA const & params,
                             ARRAY_1D_TO_CONST const & speciesConcentration,
                             ARRAY_1D & reactionRates,
                             ARRAY_2D & reactionRatesDerivatives );

  /**
   * @brief Compute the kinetic species rates for a given set of kinetic reactions.
   * @tparam PARAMS_DATA The type of the parameters data.
   * @tparam CALCULATE_DERIVATIVES Whether to calculate the derivatives of the species rates with respect to the species concentrations.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @tparam ARRAY_1D The type of the array of species rates.
   * @tparam ARRAY_2D The type of the array of species rates derivatives.
   * @param temperature The temperature of the system.
   * @param params The parameters data.
   * @param speciesConcentration The array of species concentrations.
   * @param speciesRates The array of species rates.
   * @param speciesRatesDerivatives The array of species rates derivatives.
   */
  template< typename PARAMS_DATA,
            bool CALCULATE_DERIVATIVES,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeSpeciesRates_impl( RealType const & temperature,
                            PARAMS_DATA const & params,
                            ARRAY_1D_TO_CONST const & speciesConcentration,
                            ARRAY_1D & speciesRates,
                            ARRAY_2D & speciesRatesDerivatives );

};

} // namespace bulkGeneric
} // namespace hpcReact

#include "KineticReactions_impl.hpp"
#include "common/macrosCleanup.hpp"
