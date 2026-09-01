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
#include "common/nonlinearSolvers.hpp"
#include "constitutive/activity/activity.hpp"
#include <math.h>
#include <functional>
#include <iostream>

namespace hpcReact
{
namespace massActions
{

namespace massActions_impl
{

/**
 * @brief Mass action for the secondary species, with the derivative reported through a callback.
 * @param logPrimaryActivities log of the primary species activities.
 * @param logSecondaryActivities [out] log(a_j) = -log(K_j) + nu_jw log(a_w) + sum_k nu_jk log(a_k).
 * @param derivativeFunc invoked as (j, k, nu_jk) for every secondary/primary pair, where
 *   nu_jk = d log(a_sec,j)/d log(a_prim,k). Pass a no-op to discard it.
 * @param logWaterActivity ln(a_w), picked up in proportion to the water stoichiometry. Defaults to
 *   0, i.e. an ideal solvent.
 *
 * The derivative is just the stoichiometric coefficient the loop already reads, so routing it
 * through a callback lets the value-only and value-plus-derivative public overloads share one body.
 * Callers use those overloads, not this.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename FUNC >
HPCREACT_HOST_DEVICE
inline
void calculateLogSecondaryActivities( PARAMS_DATA const & params,
                                      ARRAY_1D_TO_CONST const & logPrimaryActivities,
                                      ARRAY_1D & logSecondaryActivities,
                                      FUNC && derivativeFunc,
                                      REAL_TYPE const logWaterActivity = 0.0 )
{
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  static constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();

  for( INDEX_TYPE i = 0; i < numSecondarySpecies; ++i )
  {
    logSecondaryActivities[i] = 0.0;
  }

  for( int j=0; j<numSecondarySpecies; ++j )
  {
    // Water is not one of the species, so its stoichiometry is carried separately. The term
    // vanishes for an ideal solvent, where logWaterActivity is 0.
    logSecondaryActivities[j] = -log( params.equilibriumConstant( j ) )
                                + params.waterStoichiometry( j ) * logWaterActivity;
    for( int k=0; k<numPrimarySpecies; ++k )
    {
      logSecondaryActivities[j] += params.stoichiometricMatrix( j, k+numSecondarySpecies ) * ( logPrimaryActivities[k] );
      derivativeFunc( j, k, params.stoichiometricMatrix( j, k+numSecondarySpecies ) );
    }
  }
}

} // namespace

/**
 * @brief Secondary species activities from mass action.
 * @param logPrimaryActivities log of the primary species activities.
 * @param logSecondaryActivities [out] log(a_j) = -log(K_j) + nu_jw log(a_w) + sum_k nu_jk log(a_k).
 * @param logWaterActivity ln(a_w), picked up in proportion to the water stoichiometry. Defaults to
 *   0, i.e. an ideal solvent.
 *
 * Updates the secondary activities directly from the primary activities, with no derivative output;
 * use calculateLogSecondaryActivitiesWrtLogA when the derivative is needed. No additional activity
 * model update is needed, since mass action is stated in activities. A no-op when there are no
 * secondary species.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D >
HPCREACT_HOST_DEVICE
inline
void calculateLogSecondaryActivities( PARAMS_DATA const & params,
                                      ARRAY_1D_TO_CONST const & logPrimaryActivities,
                                      ARRAY_1D & logSecondaryActivities,
                                      REAL_TYPE const logWaterActivity = 0.0 )
{
  if constexpr( PARAMS_DATA::numSecondarySpecies() <= 0 )
  {
    return;
  }

  massActions_impl::calculateLogSecondaryActivities< REAL_TYPE,
                                                     INT_TYPE,
                                                     INDEX_TYPE >( params,
                                                                   logPrimaryActivities,
                                                                   logSecondaryActivities,
                                                                   []( INDEX_TYPE, INDEX_TYPE, REAL_TYPE ){},
                                                                   logWaterActivity );
}


/**
 * @brief Secondary species activities from mass action, and their derivatives.
 * @param logPrimaryActivities log of the primary species activities.
 * @param logSecondaryActivities [out] log(a_j) = -log(K_j) + nu_jw log(a_w) + sum_k nu_jk log(a_k).
 * @param dLogSecondaryActivities_dLogPrimaryActivities [out] d log(a_sec,j)/d log(a_prim,k), which
 *   for mass action is just the stoichiometric matrix nu.
 * @param logWaterActivity ln(a_w), picked up in proportion to the water stoichiometry. Defaults to
 *   0, i.e. an ideal solvent, and is held fixed, so it contributes nothing to the derivative.
 *
 * Same as calculateLogSecondaryActivities, but also returns the derivative.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateLogSecondaryActivitiesWrtLogA( PARAMS_DATA const & params,
                                             ARRAY_1D_TO_CONST const & logPrimaryActivities,
                                             ARRAY_1D & logSecondaryActivities,
                                             ARRAY_2D & dLogSecondaryActivities_dLogPrimaryActivities,
                                             REAL_TYPE const logWaterActivity = 0.0 )
{
  massActions_impl::calculateLogSecondaryActivities< REAL_TYPE, INT_TYPE, INDEX_TYPE >( params,
                                                                                        logPrimaryActivities,
                                                                                        logSecondaryActivities,
                                                                                        [&]( const int j, const int k, REAL_TYPE const value )
  {
    dLogSecondaryActivities_dLogPrimaryActivities[j][k] = value;
  },
                                                                                        logWaterActivity );
}

/**
 * @brief Secondary species concentrations from mass action, without updating the activity
 *   coefficients.
 * @param logPrimaryActivities log of the primary species activities.
 * @param logSecondaryActivityCoefficients log of the activity coefficients for the secondary
 *   species. Mass action yields activities; this converts them to concentrations via
 *   log(C_j) = log(a_j) - logSecondaryActivityCoefficients[j]. Pass zeros for an ideal solution,
 *   in which case the two coincide.
 * @param logWaterActivity ln(a_w), which reactions pick up in proportion to their water
 *   stoichiometry. Defaults to 0, i.e. an ideal solvent.
 *
 * The activity coefficients are taken as given, so the result is not self-consistent: the
 * concentrations returned imply an ionic strength that need not reproduce them. Use
 * calculateLogSecondarySpeciesConcentration to close that loop; this one is its ideal seed.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_1D >
HPCREACT_HOST_DEVICE
inline
void calculateLogSecondarySpeciesConcentrationNoActivityUpdate( PARAMS_DATA const & params,
                                                                ARRAY_1D_TO_CONST const & logPrimaryActivities,
                                                                ARRAY_1D_TO_CONST2 const & logSecondaryActivityCoefficients,
                                                                ARRAY_1D & logSecondarySpeciesConcentrations,
                                                                REAL_TYPE const logWaterActivity = 0.0 )
{
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();

  calculateLogSecondaryActivities< REAL_TYPE,
                                   INT_TYPE,
                                   INDEX_TYPE >( params,
                                                 logPrimaryActivities,
                                                 logSecondarySpeciesConcentrations,
                                                 logWaterActivity );

  for( INDEX_TYPE j = 0; j < numSecondarySpecies; ++j )
  {
    logSecondarySpeciesConcentrations[j] -= logSecondaryActivityCoefficients[j];
  }
}

/**
 * @copydoc calculateLogSecondarySpeciesConcentrationNoActivityUpdate
 * @param dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations [out] derivatives at
 *   fixed activity coefficients, which is the stoichiometric matrix nu. Use
 *   calculateLogSecondarySpeciesConcentrationWrtLogC when the activity coefficients vary with
 *   concentration.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_1D,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateLogSecondarySpeciesConcentrationWrtLogCNoActivityUpdate( PARAMS_DATA const & params,
                                                                       ARRAY_1D_TO_CONST const & logPrimaryActivities,
                                                                       ARRAY_1D_TO_CONST2 const & logSecondaryActivityCoefficients,
                                                                       ARRAY_1D & logSecondarySpeciesConcentrations,
                                                                       ARRAY_2D & dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                                                                       REAL_TYPE const logWaterActivity = 0.0 )
{
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  static constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();

  REAL_TYPE dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesActivities[numSecondarySpecies][numPrimarySpecies] = {{ 0.0 }};

  calculateLogSecondaryActivitiesWrtLogA< REAL_TYPE,
                                          INT_TYPE,
                                          INDEX_TYPE >( params,
                                                        logPrimaryActivities,
                                                        logSecondarySpeciesConcentrations,
                                                        dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesActivities,
                                                        logWaterActivity );

  // The activity coefficients are frozen here, so log(a) = log(C) + const and the two derivatives
  // coincide. If they vary with concentration, use calculateLogSecondarySpeciesConcentrationWrtLogC.
  for( INDEX_TYPE j = 0; j < numSecondarySpecies; ++j )
  {
    logSecondarySpeciesConcentrations[j] -= logSecondaryActivityCoefficients[j];
    for( INDEX_TYPE k = 0; k < numPrimarySpecies; ++k )
    {
      dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[j][k] =
        dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesActivities[j][k];
    }
  }
}

/**
 * @brief Secondary species concentrations, solved self-consistently with the activity model.
 * @tparam ACTIVITY_MODEL the activity model. Mass action is stated in activities, so it is
 *   activity-model dependent.
 * @param logPrimarySpeciesConcentrations log of the primary species concentrations, held fixed.
 * @param logSecondarySpeciesConcentrations [out] log of the secondary species concentrations.
 * @param logActivityCoefficients [out] log of the activity coefficients for all species, secondary
 *   rows first.
 * @param logActivities [out] log of the activities for all species.
 * @param dLogActivities_dLogSpeciesConcentrations [out] d log(a_i)/d log(c_j).
 * @param dLogActivityCoefficients_dLogSpeciesConcentrations [out] d log(activityCoefficient_i)/d log(c_j).
 * @param logWaterActivity [out] ln(a_w) at the converged state.
 * @param dLogWaterActivity_dLogSpeciesConcentrations [out] d ln(a_w)/d log(c_j).
 * @return whether the solve converged.
 *
 * The activity coefficients depend on ionic strength, which depends on the secondary concentrations
 * this function produces. An inner Newton solve on log(C_sec) closes that loop, so on return the
 * concentrations, activity coefficients, activities and a_w are mutually consistent at the given
 * primary concentrations. Note the *primary* activity coefficients move too, since ionic strength
 * does.
 *
 * Requires numSecondarySpecies > 0; the caller decides whether there is anything to solve.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_SECONDARY,
          typename ARRAY_1D_GAMMA,
          typename ARRAY_1D_ACTIVITIES,
          typename ARRAY_2D_ACTIVITIES,
          typename ARRAY_2D_GAMMA,
          typename ARRAY_1D_WATER >
HPCREACT_HOST_DEVICE
inline
bool calculateLogSecondarySpeciesConcentration( PARAMS_DATA const & params,
                                                typename ACTIVITY_MODEL::Params const & activityParams,
                                                ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                ARRAY_1D_SECONDARY & logSecondarySpeciesConcentrations,
                                                ARRAY_1D_GAMMA & logActivityCoefficients,
                                                ARRAY_1D_ACTIVITIES & logActivities,
                                                ARRAY_2D_ACTIVITIES & dLogActivities_dLogSpeciesConcentrations,
                                                ARRAY_2D_GAMMA & dLogActivityCoefficients_dLogSpeciesConcentrations,
                                                REAL_TYPE & logWaterActivity,
                                                ARRAY_1D_WATER & dLogWaterActivity_dLogSpeciesConcentrations )
{
  static constexpr int numSpecies          = PARAMS_DATA::numSpecies();
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  static constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();

  static_assert( numSecondarySpecies > 0,
                 "no secondary species to solve for; guard the call site" );
  static_assert( ACTIVITY_MODEL::Params::numSpecies() == numSpecies,
                 "activity model and reaction parameters disagree on the species count" );
  static_assert( LOGE_CONCENTRATION,
                 "only LOGE_CONCENTRATION == true is available; LOGE_CONCENTRATION == false will be "
                 "implemented upon request." );

  REAL_TYPE logSpeciesConcentration[numSpecies] = { 0.0 };
  for( INDEX_TYPE k = 0; k < numPrimarySpecies; ++k )
  {
    logSpeciesConcentration[k + numSecondarySpecies] = logPrimarySpeciesConcentrations[k];
  }

  // Initial guess from an ideal-solution pass, where the primary concentrations are also the
  // primary activities. Already exact when ACTIVITY_MODEL is Identity, so the solve below then
  // converges on its first residual evaluation.
  REAL_TYPE const logSecondaryActivityCoefficientsInitialGuess[numSecondarySpecies] = { 0.0 };
  REAL_TYPE logSecondarySpeciesConcentrationsSolution[numSecondarySpecies] = { 0.0 };
  calculateLogSecondarySpeciesConcentrationNoActivityUpdate< REAL_TYPE,
                                                             INT_TYPE,
                                                             INDEX_TYPE >( params,
                                                                           logPrimarySpeciesConcentrations,
                                                                           logSecondaryActivityCoefficientsInitialGuess,
                                                                           logSecondarySpeciesConcentrationsSolution );

  auto residualAndJacobian = [&] ( REAL_TYPE const (&x)[numSecondarySpecies],
                                   REAL_TYPE (& residual)[numSecondarySpecies],
                                   REAL_TYPE (& jacobian)[numSecondarySpecies][numSecondarySpecies] )
  {
    for( INDEX_TYPE j = 0; j < numSecondarySpecies; ++j )
    {
      logSpeciesConcentration[j] = x[j];
    }

    calculateActivities< REAL_TYPE,
                         INT_TYPE,
                         INDEX_TYPE,
                         ACTIVITY_MODEL,
                         LOGE_CONCENTRATION >( activityParams,
                                               logSpeciesConcentration,
                                               logActivities,
                                               dLogActivities_dLogSpeciesConcentrations,
                                               logActivityCoefficients,
                                               dLogActivityCoefficients_dLogSpeciesConcentrations,
                                               logWaterActivity,
                                               dLogWaterActivity_dLogSpeciesConcentrations );

    // Equilibrium constraint: log(a_j) + log(K_j) - sum_k nu_jk log(a_k) - nu_jw log(a_w) = 0.
    for( INDEX_TYPE j = 0; j < numSecondarySpecies; ++j )
    {
      REAL_TYPE const nu_jw = params.waterStoichiometry( j );

      REAL_TYPE r = logActivities[j] + log( params.equilibriumConstant( j ) )
                    - nu_jw * logWaterActivity;
      for( INDEX_TYPE k = 0; k < numPrimarySpecies; ++k )
      {
        r -= params.stoichiometricMatrix( j, k + numSecondarySpecies ) *
             logActivities[k + numSecondarySpecies];
      }
      residual[j] = r;

      // d residual_j / d log(C_sec,m). The unknown enters through the activity coefficients and
      // through a_w.
      for( INDEX_TYPE m = 0; m < numSecondarySpecies; ++m )
      {
        REAL_TYPE value = dLogActivityCoefficients_dLogSpeciesConcentrations[j][m]
                          - nu_jw * dLogWaterActivity_dLogSpeciesConcentrations[m];
        for( INDEX_TYPE k = 0; k < numPrimarySpecies; ++k )
        {
          value -= params.stoichiometricMatrix( j, k + numSecondarySpecies ) *
                   dLogActivityCoefficients_dLogSpeciesConcentrations[k + numSecondarySpecies][m];
        }
        jacobian[j][m] = ( j == m ? 1.0 : 0.0 ) + value;
      }
    }
  };

  bool const isConverged =
    nonlinearSolvers::newtonRaphson< numSecondarySpecies >( logSecondarySpeciesConcentrationsSolution,
                                                            residualAndJacobian );

  // Report the last iterate of each secondary species at nonconvergence.
  // LCOV_EXCL_START
  if( !isConverged )
  {
    printf( "calculateLogSecondarySpeciesConcentration: no convergence\n" );
    for( INDEX_TYPE j = 0; j < numSecondarySpecies; ++j )
    {
      printf( "  secondary species %2d: log c = %16.10g\n",
              static_cast< int >( j ),
              static_cast< double >( logSecondarySpeciesConcentrationsSolution[j] ) );
    }
  }
  // LCOV_EXCL_STOP

  for( INDEX_TYPE j = 0; j < numSecondarySpecies; ++j )
  {
    logSecondarySpeciesConcentrations[j] = logSecondarySpeciesConcentrationsSolution[j];
  }
  return isConverged;
}

/**
 * @copydoc calculateLogSecondarySpeciesConcentration
 * @param dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations [out]
 *   d log(C_sec)/d log(C_prim) at the converged state.
 *
 * Differentiating the converged equilibrium constraint with respect to log(C_prim) gives
 * ( I - A ) X = B, where X is the derivative to solve and I - A is the same matrix the inner
 * Newton used as its Jacobian:
 *
 *   X_jn = d log(C_sec,j)/d log(C_prim,n)                 (S x P)
 *   A_jm = sum_k nu_jk G[k+S][m]      - G[j][m]           (S x S)
 *   B_jn = nu_jn + sum_k nu_jk G[k+S][n+S] - G[j][n+S]    (S x P)
 *
 * with G = d log(activityCoefficient)/d log(C) over all species, secondary rows first.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_SECONDARY,
          typename ARRAY_1D_GAMMA,
          typename ARRAY_1D_ACTIVITIES,
          typename ARRAY_2D_ACTIVITIES,
          typename ARRAY_2D_GAMMA,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
bool calculateLogSecondarySpeciesConcentrationWrtLogC( PARAMS_DATA const & params,
                                                       typename ACTIVITY_MODEL::Params const & activityParams,
                                                       ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                       ARRAY_1D_SECONDARY & logSecondarySpeciesConcentrations,
                                                       ARRAY_1D_GAMMA & logActivityCoefficients,
                                                       ARRAY_1D_ACTIVITIES & logActivities,
                                                       ARRAY_2D_ACTIVITIES & dLogActivities_dLogSpeciesConcentrations,
                                                       ARRAY_2D_GAMMA & dLogActivityCoefficients_dLogSpeciesConcentrations,
                                                       ARRAY_2D & dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numSpecies          = PARAMS_DATA::numSpecies();
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  static constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();

  static_assert( LOGE_CONCENTRATION,
                 "only LOGE_CONCENTRATION == true is available; LOGE_CONCENTRATION == false will be "
                 "implemented upon request." );

  REAL_TYPE logWaterActivity;
  REAL_TYPE dLogWaterActivity_dLogSpeciesConcentrations[numSpecies] = { 0.0 };

  bool const isConverged =
    calculateLogSecondarySpeciesConcentration< REAL_TYPE,
                                               INT_TYPE,
                                               INDEX_TYPE,
                                               ACTIVITY_MODEL,
                                               LOGE_CONCENTRATION >( params,
                                                                     activityParams,
                                                                     logPrimarySpeciesConcentrations,
                                                                     logSecondarySpeciesConcentrations,
                                                                     logActivityCoefficients,
                                                                     logActivities,
                                                                     dLogActivities_dLogSpeciesConcentrations,
                                                                     dLogActivityCoefficients_dLogSpeciesConcentrations,
                                                                     logWaterActivity,
                                                                     dLogWaterActivity_dLogSpeciesConcentrations );

  REAL_TYPE identityMinusA[numSecondarySpecies][numSecondarySpecies] = {{ 0.0 }};
  REAL_TYPE rhs[numSecondarySpecies][numPrimarySpecies] = {{ 0.0 }};

  for( INDEX_TYPE j = 0; j < numSecondarySpecies; ++j )
  {
    REAL_TYPE const nu_jw = params.waterStoichiometry( j );

    for( INDEX_TYPE m = 0; m < numSecondarySpecies; ++m )
    {
      REAL_TYPE a = -dLogActivityCoefficients_dLogSpeciesConcentrations[j][m]
                    + nu_jw * dLogWaterActivity_dLogSpeciesConcentrations[m];
      for( INDEX_TYPE k = 0; k < numPrimarySpecies; ++k )
      {
        a += params.stoichiometricMatrix( j, k + numSecondarySpecies ) *
             dLogActivityCoefficients_dLogSpeciesConcentrations[k + numSecondarySpecies][m];
      }
      identityMinusA[j][m] = ( j == m ? 1.0 : 0.0 ) - a;
    }

    for( INDEX_TYPE n = 0; n < numPrimarySpecies; ++n )
    {
      REAL_TYPE b = params.stoichiometricMatrix( j, n + numSecondarySpecies ) -
                    dLogActivityCoefficients_dLogSpeciesConcentrations[j][n + numSecondarySpecies]
                    + nu_jw * dLogWaterActivity_dLogSpeciesConcentrations[n + numSecondarySpecies];
      for( INDEX_TYPE k = 0; k < numPrimarySpecies; ++k )
      {
        b += params.stoichiometricMatrix( j, k + numSecondarySpecies ) *
             dLogActivityCoefficients_dLogSpeciesConcentrations[k + numSecondarySpecies][n + numSecondarySpecies];
      }
      rhs[j][n] = b;
    }
  }

  for( INDEX_TYPE n = 0; n < numPrimarySpecies; ++n )
  {
    REAL_TYPE matrix[numSecondarySpecies][numSecondarySpecies];
    REAL_TYPE column[numSecondarySpecies];
    REAL_TYPE solution[numSecondarySpecies];

    for( INDEX_TYPE j = 0; j < numSecondarySpecies; ++j )
    {
      column[j] = rhs[j][n];
      for( INDEX_TYPE m = 0; m < numSecondarySpecies; ++m )
      {
        matrix[j][m] = identityMinusA[j][m];
      }
    }

    solveNxN_pivoted< REAL_TYPE, numSecondarySpecies >( matrix, column, solution );

    for( INDEX_TYPE j = 0; j < numSecondarySpecies; ++j )
    {
      dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[j][n] = solution[j];
    }
  }

  return isConverged;
}

/**
 * @brief Aggregate (total) primary concentrations and their derivatives, in log-concentration space.
 * @param logSecondarySpeciesConcentrations log of the secondary concentrations, already solved
 *   consistently with the activity model.
 * @param dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations d log(C_sec)/d log(C_prim)
 *   for that same converged state.
 *
 * A pure mole balance, with i and k indexing primary species, j secondary:
 *   T_i                  = C_prim,i + sum_j nu_ji C_sec,j
 *   dT_i/dlog(C_prim,k)  = delta_ik C_prim,i + sum_j nu_ji C_sec,j X_jk
 * where X_jk = d log(C_sec,j)/d log(C_prim,k).
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D_TO_CONST,
          typename ARRAY_1D_PRIMARY,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateAggregatePrimaryConcentrationsWrtLogC( PARAMS_DATA const & params,
                                                     ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                     ARRAY_1D_TO_CONST2 const & logSecondarySpeciesConcentrations,
                                                     ARRAY_2D_TO_CONST const & dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                                                     ARRAY_1D_PRIMARY & aggregatePrimarySpeciesConcentrations,
                                                     ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();

  for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
  {
    for( INDEX_TYPE j = 0; j < numPrimarySpecies; ++j )
    {
      dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations[i][j] = 0.0;
    }
  }

  for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
  {
    REAL_TYPE const primarySpeciesConcentration_i = exp( logPrimarySpeciesConcentrations[i] );
    aggregatePrimarySpeciesConcentrations[i] = primarySpeciesConcentration_i;
    dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, i ) = primarySpeciesConcentration_i;
  }

  if constexpr( numSecondarySpecies > 0 )
  {
    for( INDEX_TYPE j = 0; j < numSecondarySpecies; ++j )
    {
      REAL_TYPE const secondarySpeciesConcentration_j = exp( logSecondarySpeciesConcentrations[j] );

      for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
      {
        REAL_TYPE const nu_ji = params.stoichiometricMatrix( j, i + numSecondarySpecies );

        aggregatePrimarySpeciesConcentrations[i] += nu_ji * secondarySpeciesConcentration_j;

        for( INDEX_TYPE k = 0; k < numPrimarySpecies; ++k )
        {
          dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, k ) +=
            nu_ji * secondarySpeciesConcentration_j *
            dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[j][k];
        }
      }
    }
  }
  else
  {
    HPCREACT_UNUSED_VAR( params );
    HPCREACT_UNUSED_VAR( logSecondarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations );
  }
}

/**
 * @copydoc calculateAggregatePrimaryConcentrationsWrtLogC
 *
 * Also accumulates the mobile-only aggregate, excluding immobile secondary species.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename PARAMS_DATA,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D_TO_CONST,
          typename ARRAY_1D_PRIMARY,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE
inline
void calculateTotalAndMobileAggregatePrimaryConcentrationsWrtLogC( PARAMS_DATA const & params,
                                                                   ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentrations,
                                                                   ARRAY_1D_TO_CONST2 const & logSecondarySpeciesConcentrations,
                                                                   ARRAY_2D_TO_CONST const & dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                                                                   ARRAY_1D_PRIMARY & aggregatePrimarySpeciesConcentrations,
                                                                   ARRAY_1D_PRIMARY & mobileAggregatePrimarySpeciesConcentrations,
                                                                   ARRAY_2D & dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations,
                                                                   ARRAY_2D & dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations )
{
  static constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();

  for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
  {
    for( INDEX_TYPE j = 0; j < numPrimarySpecies; ++j )
    {
      dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations[i][j] = 0.0;
      dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations[i][j] = 0.0;
    }
  }

  for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
  {
    REAL_TYPE const primarySpeciesConcentration_i = exp( logPrimarySpeciesConcentrations[i] );
    aggregatePrimarySpeciesConcentrations[i] = primarySpeciesConcentration_i;
    mobileAggregatePrimarySpeciesConcentrations[i] = primarySpeciesConcentration_i;
    dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, i ) = primarySpeciesConcentration_i;
    dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, i ) = primarySpeciesConcentration_i;
  }

  if constexpr( numSecondarySpecies > 0 )
  {
    for( INDEX_TYPE j = 0; j < numSecondarySpecies; ++j )
    {
      REAL_TYPE const secondarySpeciesConcentration_j = exp( logSecondarySpeciesConcentrations[j] );
      REAL_TYPE const mobileFlag = params.mobileSecondarySpeciesFlag( j );

      for( INDEX_TYPE i = 0; i < numPrimarySpecies; ++i )
      {
        REAL_TYPE const nu_ji = params.stoichiometricMatrix( j, i + numSecondarySpecies );

        aggregatePrimarySpeciesConcentrations[i] += nu_ji * secondarySpeciesConcentration_j;
        mobileAggregatePrimarySpeciesConcentrations[i] += nu_ji * secondarySpeciesConcentration_j * mobileFlag;

        for( INDEX_TYPE k = 0; k < numPrimarySpecies; ++k )
        {
          // What secondary species j adds to the aggregate derivative for primary species i.
          REAL_TYPE const aggregateDerivativeContribution =
            nu_ji * secondarySpeciesConcentration_j *
            dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[j][k];

          dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, k ) +=
            aggregateDerivativeContribution;
          dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, k ) +=
            aggregateDerivativeContribution * mobileFlag;
        }
      }
    }
  }
  else
  {
    HPCREACT_UNUSED_VAR( params );
    HPCREACT_UNUSED_VAR( logSecondarySpeciesConcentrations );
    HPCREACT_UNUSED_VAR( dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations );
  }
}


} // namespace massActions
} // namespace hpcReact
