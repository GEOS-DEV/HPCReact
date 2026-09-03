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
#include "common/CArrayWrapper.hpp"
#include "common/DirectSystemSolve.hpp"
#include "common/printers.hpp"
#include "constitutive/activity/activity.hpp"

#include <iostream>


namespace hpcReact
{
namespace reactionsSystems
{

/**
 * @brief This class implements components required to calculate equilibrium
 *        reactions for a given set of species. The class also proovides
 *        device callable methods to enforce equilibrium pointwise and also
 *        provides methods to launch batch processing of material points.
 * @tparam REAL_TYPE The type of the real numbers used in the class.
 * @tparam INT_TYPE The type of the integers used in the class.
 * @tparam INDEX_TYPE The type of the indices used in the class.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL >
class EquilibriumReactions
{
public:
  /// alias for type of the real numbers used in the class.
  using RealType = REAL_TYPE;

  /// alias for type of the integers used in the class.
  using IntType = INT_TYPE;

  /// alias for type of the indices used in the class.
  using IndexType = INDEX_TYPE;



  /**
   * @brief This method enforces equilibrium for a given set of species using
   *        reaction extents.
   * @param temperature The temperature of the system.
   * @param params The parameters for the equilibrium reactions.
   * @param activityParams The parameters for the activity model.
   * @param speciesConcentration0 The initial species concentrations.
   * @param speciesConcentration The species concentrations to be updated.
   * @details This method uses the reaction extents to enforce equilibrium
   *          for a given set of species. It uses the computeResidualAndJacobian
   *          method to compute the residual and jacobian for the system and
   *          then uses a direct solver to solve the system. The solution is
   *          then used to update the species concentrations.
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST >
  static HPCREACT_HOST_DEVICE
  void
  enforceEquilibrium_Extents( RealType const & temperature,
                              PARAMS_DATA const & params,
                              typename ACTIVITY_MODEL::Params const & activityParams,
                              ARRAY_1D_TO_CONST const & speciesConcentration0,
                              ARRAY_1D & speciesConcentration );

  /**
   * @brief This method enforces equilibrium for a given set of species using
   *        aggregate primary concentrations.
   * @param temperature The temperature of the system.
   * @param params The parameters for the equilibrium reactions.
   * @param activityParams The parameters for the activity model.
   * @param speciesConcentration0 The initial species concentrations.
   * @param speciesConcentration The species concentrations to be updated.
   * @return whether the solve converged.
   * @details This method uses the aggregate primary concentrations to enforce
   *          equilibrium for a given set of species. It uses the
   *          computeResidualAndJacobianAggregatePrimaryConcentrations method to
   *          compute the residual and jacobian for the system and then uses a
   *          direct solver to solve the system. The solution is then used to
   *          update the species concentrations.
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST >
  static HPCREACT_HOST_DEVICE
  bool
  enforceEquilibrium_LogAggregate( RealType const & temperature,
                                   PARAMS_DATA const & params,
                                   typename ACTIVITY_MODEL::Params const & activityParams,
                                   ARRAY_1D_TO_CONST const & speciesConcentration0,
                                   ARRAY_1D & speciesConcentration );


  /**
   * @brief This method enforces equilibrium for a given set of species using
   *        log of aggregate primary concentrations.
   * @tparam PARAMS_DATA The type of the parameters data.
   * @tparam ARRAY_1D The type of the array of species concentrations.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @param temperature The temperature of the system.
   * @param params The parameters for the equilibrium reactions.
   * @param activityParams The parameters for the activity model.
   * @param targetAggregatePrimarySpeciesConcentration The target aggregate
   *        primary species concentration.
   * @param logPrimarySpeciesConcentration0 The initial value of the log of
   *        the primary species concentrations.
   * @param logPrimarySpeciesConcentration [out] The log of the primary species concentrations.
   * @param logSecondarySpeciesConcentration [out] The log of the secondary species concentrations
   *        at the converged state, which the last residual evaluation produces anyway.
   * @return whether the solve converged, both the outer Newton loop and every speciation solve
   *         inside it.
   * @details This method uses the log of aggregate primary concentrations to enforce
   *          equilibrium for a given set of species. It uses the
   *          computeResidualAndJacobianLogAggregate method to compute the residual and
   *          jacobian for the system and then uses a direct solver to solve the system.
   *          The solution is then used to update the species concentrations.
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_SECONDARY >
  static HPCREACT_HOST_DEVICE
  bool
  enforceEquilibrium_Aggregate( RealType const & temperature,
                                PARAMS_DATA const & params,
                                typename ACTIVITY_MODEL::Params const & activityParams,
                                ARRAY_1D_TO_CONST const & targetAggregatePrimarySpeciesConcentration,
                                ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration0,
                                ARRAY_1D & logPrimarySpeciesConcentration,
                                ARRAY_1D_SECONDARY & logSecondarySpeciesConcentration );

  /**
   * @brief Overload for callers that want only the primary concentrations.
   * @tparam PARAMS_DATA The type of the parameters data.
   * @tparam ARRAY_1D The type of the array of species concentrations.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @param temperature The temperature of the system.
   * @param params The parameters for the equilibrium reactions.
   * @param activityParams The parameters for the activity model.
   * @param targetAggregatePrimarySpeciesConcentration The target aggregate
   *        primary species concentration.
   * @param logPrimarySpeciesConcentration0 The initial value of the log of
   *        the primary species concentrations.
   * @param logPrimarySpeciesConcentration [out] The log of the primary species concentrations.
   * @return whether the solve converged.
   * @details The secondary concentrations are computed either way, so prefer the form above over
   *          recovering them with a second speciation solve.
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST >
  static HPCREACT_HOST_DEVICE
  bool
  enforceEquilibrium_Aggregate( RealType const & temperature,
                                PARAMS_DATA const & params,
                                typename ACTIVITY_MODEL::Params const & activityParams,
                                ARRAY_1D_TO_CONST const & targetAggregatePrimarySpeciesConcentration,
                                ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration0,
                                ARRAY_1D & logPrimarySpeciesConcentration )
  {
    static constexpr INDEX_TYPE numSecondarySpeciesStorage =
      PARAMS_DATA::numSecondarySpecies() > 0 ? PARAMS_DATA::numSecondarySpecies() : 1;

    RealType logSecondarySpeciesConcentration[numSecondarySpeciesStorage] = { 0.0 };

    return enforceEquilibrium_Aggregate( temperature,
                                         params,
                                         activityParams,
                                         targetAggregatePrimarySpeciesConcentration,
                                         logPrimarySpeciesConcentration0,
                                         logPrimarySpeciesConcentration,
                                         logSecondarySpeciesConcentration );
  }

  /**
   * @brief This method computes the residual and jacobian when using reaction extents to solve
   *       for the equilibrium of a given set of species.
   * @tparam PARAMS_DATA The type of the parameters data.
   * @tparam ARRAY_1D The type of the array of species concentrations.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @tparam ARRAY_1D_TO_CONST2 The type of the array of reaction extents.
   * @tparam ARRAY_2D The type of the array of jacobian.
   * @param temperature The temperature of the system.
   * @param params The parameters for the equilibrium reactions.
   * @param activityParams The parameters for the activity model.
   * @param speciesConcentration0 The initial species concentrations.
   * @param xi The reaction extents.
   * @param residual The residual.
   * @param jacobian The jacobian.
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST2,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeResidualAndJacobianReactionExtents( RealType const & temperature,
                                             PARAMS_DATA const & params,
                                             typename ACTIVITY_MODEL::Params const & activityParams,
                                             ARRAY_1D_TO_CONST const & speciesConcentration0,
                                             ARRAY_1D_TO_CONST2 const & xi,
                                             ARRAY_1D & residual,
                                             ARRAY_2D & jacobian );

  /**
   * @brief This method computes the residual and jacobian when using aggregate primary concentrations
   *        to solve for the equilibrium of a given set of species.
   * @tparam PARAMS_DATA The type of the parameters data.
   * @tparam ARRAY_1D The type of the array of species concentrations.
   * @tparam ARRAY_1D_TO_CONST The type of the array of target aggregate primary concentrations.
   * @tparam ARRAY_1D_TO_CONST The type of the array of log primary species concentrations.
   * @tparam ARRAY_2D The type of the array of jacobian.
   * @param temperature The temperature of the system.
   * @param params The parameters for the equilibrium reactions.
   * @param activityParams The parameters for the activity model.
   * @param targetAggregatePrimaryConcentrations The target aggregate primary concentrations.
   * @param logPrimarySpeciesConcentration The log of the primary species concentrations.
   * @param residual The residual.
   * @param jacobian The jacobian.
   * @param logSecondarySpeciesConcentration [out] The log of the secondary species concentrations
   *        the speciation solve produced at these primary concentrations.
   * @return whether the inner speciation solve converged.
   */
  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST2,
            typename ARRAY_2D,
            typename ARRAY_1D_SECONDARY >
  static HPCREACT_HOST_DEVICE bool
  computeResidualAndJacobianAggregatePrimaryConcentrations( RealType const & temperature,
                                                            PARAMS_DATA const & params,
                                                            typename ACTIVITY_MODEL::Params const & activityParams,
                                                            ARRAY_1D_TO_CONST const & targetAggregatePrimaryConcentrations,
                                                            ARRAY_1D_TO_CONST2 const & logPrimarySpeciesConcentration,
                                                            ARRAY_1D & residual,
                                                            ARRAY_2D & jacobian,
                                                            ARRAY_1D_SECONDARY & logSecondarySpeciesConcentration );
};



} // namespace reactionsSystems
} // namespace hpcReact

#if !defined(__INTELLISENSE__)
#include "EquilibriumReactionsAggregatePrimaryConcentration_impl.hpp"
#include "EquilibriumReactionsReactionExtents_impl.hpp"
#endif
