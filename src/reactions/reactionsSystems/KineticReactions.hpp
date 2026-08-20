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
#include "constitutive/activity/activity.hpp"
#include "Parameters.hpp"

#include <stdexcept>

/** @file KineticReactions.hpp
 *  @brief Header file for the KineticReactions class.
 *  @author HPC-REACT Team
 *  @date 2023
 */

namespace hpcReact
{
namespace reactionsSystems
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
          typename ACTIVITY_MODEL,
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
   * @copydoc KineticReactions::computeReactionRates()
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        typename ACTIVITY_MODEL::Params const & activityParams,
                        ARRAY_1D_TO_CONST const & speciesConcentration,
                        ARRAY_1D & reactionRates,
                        ARRAY_2D & dReactionRates_dConcentration )
  {

    RealType activities[ PARAMS_DATA::numSpecies() ];
    RealType dActivities_dConcentration[ PARAMS_DATA::numSpecies() ][ PARAMS_DATA::numSpecies() ] = { 0.0 };
    RealType dReactionRates_dActivities[ PARAMS_DATA::numReactions() ][ PARAMS_DATA::numSpecies() ] = { 0.0 };

    calculateActivities< RealType,
                         IntType,
                         IndexType,
                         ACTIVITY_MODEL,
                         LOGE_CONCENTRATION >( activityParams,
                                               speciesConcentration,
                                               activities,
                                               dActivities_dConcentration );

    computeReactionRatesElementary_impl< PARAMS_DATA, true >( temperature,
                                                              params,
                                                              activities,
                                                              reactionRates,
                                                              dReactionRates_dActivities );

    // chain rule to get dReactionRate_dConcentration
    for( IntType r=0; r<PARAMS_DATA::numReactions(); ++r )
    {
      for( IntType i=0; i<PARAMS_DATA::numSpecies(); ++i )
      {
        dReactionRates_dConcentration( r, i ) = 0.0;
        for( IntType j=0; j<PARAMS_DATA::numSpecies(); ++j )
        {
          dReactionRates_dConcentration( r, i ) += dReactionRates_dActivities[r][j] * dActivities_dConcentration[j][i];
        }
      }
    }
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
                        typename ACTIVITY_MODEL::Params const & activityParams,
                        ARRAY_1D_TO_CONST const & speciesConcentration,
                        ARRAY_1D & reactionRates )
  {


    RealType activities[ PARAMS_DATA::numSpecies() ];
    RealType dActivities_dConcentration[ PARAMS_DATA::numSpecies() ][ PARAMS_DATA::numSpecies() ] = { 0.0 };
    RealType dReactionRates_dActivities[ PARAMS_DATA::numReactions() ][ PARAMS_DATA::numSpecies() ] = { 0.0 };

    calculateActivities< RealType,
                         IntType,
                         IndexType,
                         ACTIVITY_MODEL,
                         LOGE_CONCENTRATION >( activityParams,
                                               speciesConcentration,
                                               activities,
                                               dActivities_dConcentration );

    computeReactionRatesElementary_impl< PARAMS_DATA, false >( temperature,
                                                               params,
                                                               activities,
                                                               reactionRates,
                                                               dReactionRates_dActivities );
    HPCREACT_UNUSED_VAR( dReactionRates_dActivities );
  }

  /**
   * @brief Compute the reaction rates for a given set of species concentrations and surface area.
   * @tparam PARAMS_DATA The type of the parameters data.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations
   * @tparam ARRAY_1D_SA The type of the array of surface area.
   * @tparam ARRAY_1D The type of the array of reaction rates.
   * @tparam ARRAY_2D The type of the array of reaction rates derivatives.
   * @param temperature The temperature of the system.
   * @param params The parameters data.
   * @param speciesConcentration The array of species concentrations.
   * @param surfaceArea The array of surface area.
   * @param reactionRates The array of reaction rates.
   * @param dReactionRates_dConcentration The array of reaction rates derivatives.
   * @details
   *   This function computes the reaction rates for a given set of reactions,
   *   taking into account the surface area of the reactions. If
   *   CALCULATE_DERIVATIVES is true, it also computes the derivatives of the
   *   reaction rates with respect to the species concentrations.
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_SA,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        typename ACTIVITY_MODEL::Params const & activityParams,
                        ARRAY_1D_TO_CONST const & speciesConcentration,
                        ARRAY_1D_SA const & surfaceArea,
                        ARRAY_1D & reactionRates,
                        ARRAY_2D & dReactionRates_dConcentration )
  {

    RealType activities[ PARAMS_DATA::numSpecies() ];
    RealType dActivities_dConcentration[ PARAMS_DATA::numSpecies() ][ PARAMS_DATA::numSpecies() ]{  };
    RealType dReactionRates_dActivities[ PARAMS_DATA::numReactions() ][ PARAMS_DATA::numSpecies() ]{  };

    calculateActivities< RealType,
                         IntType,
                         IndexType,
                         ACTIVITY_MODEL,
                         LOGE_CONCENTRATION >( activityParams,
                                               speciesConcentration,
                                               activities,
                                               dActivities_dConcentration );

    if( params.reactionRateLawOption() == ReactionRateLawOption::Elementary )
    {
      computeReactionRatesElementary_impl< PARAMS_DATA, true >( temperature,
                                                                params,
                                                                activities,
                                                                reactionRates,
                                                                dReactionRates_dActivities );
    }
    else if( params.reactionRateLawOption() == ReactionRateLawOption::Affinity )
    {
      computeReactionRatesAffinity_impl< PARAMS_DATA, true >( temperature,
                                                              params,
                                                              activities,
                                                              surfaceArea,
                                                              reactionRates,
                                                              dReactionRates_dActivities );
    }

    // chain rule to get dReactionRate_dConcentration
    for( IntType r=0; r<PARAMS_DATA::numReactions(); ++r )
    {
      for( IntType i=0; i<PARAMS_DATA::numSpecies(); ++i )
      {
        dReactionRates_dConcentration( r, i ) = 0.0;
        for( IntType j=0; j<PARAMS_DATA::numSpecies(); ++j )
        {
          dReactionRates_dConcentration( r, i ) += dReactionRates_dActivities[r][j] * dActivities_dConcentration[j][i];
        }
      }
    }
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
                       typename ACTIVITY_MODEL::Params const & activityParams,
                       ARRAY_1D_TO_CONST const & speciesConcentration,
                       ARRAY_1D & speciesRates,
                       ARRAY_2D & dSpeciesRates_dConcentration )
  {

    RealType activities[ PARAMS_DATA::numSpecies() ];
    RealType dActivities_dConcentration[ PARAMS_DATA::numSpecies() ][ PARAMS_DATA::numSpecies() ] {};
    RealType dSpeciesRates_dActivities[ PARAMS_DATA::numSpecies() ][ PARAMS_DATA::numSpecies() ] {};

    calculateActivities< RealType,
                         IntType,
                         IndexType,
                         ACTIVITY_MODEL,
                         LOGE_CONCENTRATION >( activityParams,
                                               speciesConcentration,
                                               activities,
                                               dActivities_dConcentration );

    computeSpeciesRates_impl< PARAMS_DATA, true >( temperature,
                                                   params,
                                                   activities,
                                                   speciesRates,
                                                   dSpeciesRates_dActivities );
    // chain rule to get dSpeciesRates_dConcentration
    for( IntType i=0; i<PARAMS_DATA::numSpecies(); ++i )
    {
      for( IntType j=0; j<PARAMS_DATA::numSpecies(); ++j )
      {
        dSpeciesRates_dConcentration( i, j ) = 0.0;
        for( IntType k=0; k<PARAMS_DATA::numSpecies(); ++k )
        {
          dSpeciesRates_dConcentration( i, j ) += dSpeciesRates_dActivities[i][k] * dActivities_dConcentration[k][j];
        }
      }
    }
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
                       typename ACTIVITY_MODEL::Params const & activityParams,
                       ARRAY_1D_TO_CONST const & speciesConcentration,
                       ARRAY_1D & speciesRates )
  {
    RealType activities[ PARAMS_DATA::numSpecies() ];
    RealType dActivities_dConcentration[ PARAMS_DATA::numSpecies() ][ PARAMS_DATA::numSpecies() ] = { 0.0 };

    calculateActivities< RealType,
                         IntType,
                         IndexType,
                         ACTIVITY_MODEL,
                         LOGE_CONCENTRATION >( activityParams,
                                               speciesConcentration,
                                               activities,
                                               dActivities_dConcentration );

    char speciesRatesDerivatives;

    computeSpeciesRates_impl< PARAMS_DATA, false >( temperature,
                                                    params,
                                                    activities,
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
   * @note Currently uninstantiated, and will not compile as written: the body calls
   *   computeSpeciesRates without the activityParams argument every overload now requires. Add that
   *   parameter here and pass it through before using this.
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
  computeReactionRatesElementary_impl( RealType const & temperature,
                                       PARAMS_DATA const & params,
                                       ARRAY_1D_TO_CONST const & activities,
                                       ARRAY_1D & reactionRates,
                                       ARRAY_2D & dReactionRates_dActivities );

/**
 * @brief Compute the reaction rates from the departure of the activity quotient from equilibrium,
 *   \f$ \dot{R}_r = k_r A_r ( 1 - Q_r / K_r ) \f$.
 *
 * @tparam PARAMS_DATA
 * @tparam CALCULATE_DERIVATIVES
 * @tparam ARRAY_1D_TO_CONST
 * @tparam ARRAY_1D
 * @tparam ARRAY_2D
 * @param temperature
 * @param params
 * @param activities
 * @param surfaceArea
 * @param reactionRates
 * @param dReactionRates_dActivities
 * @return HPCREACT_HOST_DEVICE
 */
  template< typename PARAMS_DATA,
            bool CALCULATE_DERIVATIVES,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_SA,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeReactionRatesAffinity_impl( RealType const & temperature,
                                     PARAMS_DATA const & params,
                                     ARRAY_1D_TO_CONST const & activities,
                                     ARRAY_1D_SA const & surfaceArea,
                                     ARRAY_1D & reactionRates,
                                     ARRAY_2D & dReactionRates_dActivities );

  /**
   * @brief Compute the kinetic species rates for a given set of kinetic reactions.
   * @tparam PARAMS_DATA The type of the parameters data.
   * @tparam CALCULATE_DERIVATIVES Whether to calculate the derivatives of the species rates with respect to the species concentrations.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @tparam ARRAY_1D The type of the array of species rates.
   * @tparam ARRAY_2D The type of the array of species rates derivatives.
   * @param temperature The temperature of the system.
   * @param params The parameters data.
   * @param activities The array of activities.
   * @param speciesRates The array of species rates.
   * @param dReactionRates_dActivities The array of species rates derivatives.
   */
  template< typename PARAMS_DATA,
            bool CALCULATE_DERIVATIVES,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeSpeciesRates_impl( RealType const & temperature,
                            PARAMS_DATA const & params,
                            ARRAY_1D_TO_CONST const & activities,
                            ARRAY_1D & speciesRates,
                            ARRAY_2D & dReactionRates_dActivities );

};

} // namespace reactionsSystems
} // namespace hpcReact

#include "KineticReactions_impl.hpp"
#include "common/macrosCleanup.hpp"
