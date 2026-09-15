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

#if defined(__INTELLISENSE__)
#include "EquilibriumReactions.hpp"
#endif

#include "common/constants.hpp"
#include "reactions/massActions/MassActions.hpp"
#include "constitutive/activity/Identity.hpp"

#include <math.h>
#include <type_traits>

namespace hpcReact
{
namespace reactionsSystems
{

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_CONSTRAINT,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D,
          typename ARRAY_1D_SECONDARY,
          typename ARRAY_1D_AGGREGATE >
HPCREACT_HOST_DEVICE
inline
bool
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE,
                      ACTIVITY_MODEL >::computeResidualAndJacobianPrimaryConcentrations( RealType const & temperature,
                                                                                         PARAMS_DATA const & params,
                                                                                         typename ACTIVITY_MODEL::Params const & activityParams,
                                                                                         ARRAY_1D_CONSTRAINT const & constraintType,
                                                                                         ARRAY_1D_TO_CONST const & constraintValue,
                                                                                         ARRAY_1D_TO_CONST2 const & logPrimarySpeciesConcentration,
                                                                                         ARRAY_1D & residual,
                                                                                         ARRAY_2D & jacobian,
                                                                                         ARRAY_1D_SECONDARY & logSecondarySpeciesConcentration,
                                                                                         ARRAY_1D_AGGREGATE & aggregatePrimarySpeciesConcentration )
{
  HPCREACT_UNUSED_VAR( temperature );
  static constexpr int numSpecies = PARAMS_DATA::numSpecies();
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  static constexpr int numSecondarySpeciesStorage = numSecondarySpecies > 0 ? numSecondarySpecies : 1;
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  bool speciationConverged = true;
  RealType aggregatePrimaryConcentrations[numPrimarySpecies] = {0.0};
  RealType mobileAggregatePrimaryConcentrations[numPrimarySpecies] = {0.0};
  ARRAY_2D dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations = {{{0.0}}};
  ARRAY_2D dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations = {{{0.0}}};
  RealType dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[numSecondarySpeciesStorage][numPrimarySpecies] = {{0.0}};
  RealType logActivities[numSpecies] = {0.0};
  RealType dLogActivities_dLogSpeciesConcentrations[numSpecies][numSpecies] = {{0.0}};

  if constexpr( numPrimarySpecies > 0 )
  {
    RealType logActivityCoefficients[numSpecies] = {0.0};
    RealType dLogActivityCoefficients_dLogSpeciesConcentrations[numSpecies][numSpecies] = {{0.0}};

    if constexpr( numSecondarySpecies > 0 )
    {
      // Secondary concentrations, activity coefficients and activities are solved together at the
      // given primary concentrations, so on return all three are mutually consistent and satisfy
      // the mass action law, and the derivative is the exact one for that converged state rather
      // than the frozen activity coefficient approximation.
      speciationConverged =
        massActions::calculateLogSecondarySpeciesConcentrationWrtLogC< REAL_TYPE,
                                                                       INT_TYPE,
                                                                       INDEX_TYPE,
                                                                       ACTIVITY_MODEL,
                                                                       true >( params,
                                                                               activityParams,
                                                                               logPrimarySpeciesConcentration,
                                                                               logSecondarySpeciesConcentration,
                                                                               logActivityCoefficients,
                                                                               logActivities,
                                                                               dLogActivities_dLogSpeciesConcentrations,
                                                                               dLogActivityCoefficients_dLogSpeciesConcentrations,
                                                                               dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations );
    }
    else
    {
      // Nothing to speciate: the primary species are all the species there are, so the activities
      // and their derivatives come straight from the activity model. Without this a pX row would
      // read ln(a) = 0 and contribute an all-zero jacobian row.
      RealType logWaterActivity = 0.0;
      RealType dLogWaterActivity_dLogSpeciesConcentrations[numSpecies] = {0.0};

      calculateActivities< REAL_TYPE,
                           INT_TYPE,
                           INDEX_TYPE,
                           ACTIVITY_MODEL,
                           true >( activityParams,
                                   logPrimarySpeciesConcentration,
                                   logActivities,
                                   dLogActivities_dLogSpeciesConcentrations,
                                   logActivityCoefficients,
                                   dLogActivityCoefficients_dLogSpeciesConcentrations,
                                   logWaterActivity,
                                   dLogWaterActivity_dLogSpeciesConcentrations );
    }
  }

  // Pure mole balance on the solved state: the activity model is already accounted for in the
  // calculation of secondary species. ChargeBalance constraint takes the mobile aggregate.
  massActions::calculateTotalAndMobileAggregatePrimaryConcentrationsWrtLogC< REAL_TYPE, INT_TYPE, INDEX_TYPE >(
    params,
    logPrimarySpeciesConcentration,
    logSecondarySpeciesConcentration,
    dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
    aggregatePrimaryConcentrations,
    mobileAggregatePrimaryConcentrations,
    dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations,
    dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );

  for( IndexType i = 0; i < numPrimarySpecies; ++i )
  {
    if( constraintType[i] == PrimarySpeciesConstraintType::pX )
    {
      // ln(a_i) + pX_i * ln(10) = 0
      residual[i] = logActivities[i + numSecondarySpecies] + constraintValue[i] * constants::ln10;

      // d ln(a_i)/d ln(C_prim,k) at the converged speciation
      for( IndexType k = 0; k < numPrimarySpecies; ++k )
      {
        RealType value = dLogActivities_dLogSpeciesConcentrations[i + numSecondarySpecies][k + numSecondarySpecies];
        for( IndexType j = 0; j < numSecondarySpecies; ++j )
        {
          value += dLogActivities_dLogSpeciesConcentrations[i + numSecondarySpecies][j] *
                   dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[j][k];
        }
        jacobian( i, k ) = -value;
      }
    }
    else if( constraintType[i] == PrimarySpeciesConstraintType::ChargeBalance )
    {
      RealType chargeSum = 0.0;
      RealType chargeScale = 0.0;
      for( IndexType p = 0; p < numPrimarySpecies; ++p )
      {
        RealType const term = activityParams.m_speciesCharge[p + numSecondarySpecies] *
                              mobileAggregatePrimaryConcentrations[p];
        chargeSum += term;
        chargeScale += fabs( term );
      }

      // Nothing to balance: every primary species is neutral, or the activity parameters were built
      // without charges. Report failure rather than assembling a singular row.
      if( !( chargeScale > 0.0 ) )
      {
        return false;
      }

      residual[i] = chargeSum / chargeScale;

      for( IndexType j = 0; j < numPrimarySpecies; ++j )
      {
        RealType dChargeSum = 0.0;
        for( IndexType p = 0; p < numPrimarySpecies; ++p )
        {
          dChargeSum += activityParams.m_speciesCharge[p + numSecondarySpecies] *
                        dMobileAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( p, j );
        }
        jacobian( i, j ) = -dChargeSum / chargeScale;
      }
    }
    // else if( constraintType[i] == PrimarySpeciesConstraintType::MineralEquilibrium ) // TODO
    else // constraintType[i] == PrimarySpeciesConstraintType::AggregateConcentration
    {
      residual[i] = -(1.0 - aggregatePrimaryConcentrations[i] / constraintValue[i]);
      for( IndexType j = 0; j < numPrimarySpecies; ++j )
      {
        jacobian( i, j ) = -dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations[i][j] / constraintValue[i];
      }
    }
  }

  for( IndexType i = 0; i < numPrimarySpecies; ++i )
  {
    aggregatePrimarySpeciesConcentration[i] = aggregatePrimaryConcentrations[i];
  }

  return speciationConverged;
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_SECONDARY >
HPCREACT_HOST_DEVICE inline
bool
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE,
                      ACTIVITY_MODEL >::enforceEquilibrium_PrimaryConcentrations( REAL_TYPE const & temperature,
                                                                                  PARAMS_DATA const & params,
                                                                                  typename ACTIVITY_MODEL::Params const & activityParams,
                                                                                  ARRAY_1D_TO_CONST const & targetAggregatePrimarySpeciesConcentration,
                                                                                  ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration0,
                                                                                  ARRAY_1D & logPrimarySpeciesConcentration,
                                                                                  ARRAY_1D_SECONDARY & logSecondarySpeciesConcentration )
{
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  PrimarySpeciesConstraintType constraintType[numPrimarySpecies];
  for( int i=0; i<numPrimarySpecies; ++i )
  {
    constraintType[i] = PrimarySpeciesConstraintType::AggregateConcentration;
  }

  return enforceEquilibrium_PrimaryConcentrations( temperature,
                                                   params,
                                                   activityParams,
                                                   constraintType,
                                                   targetAggregatePrimarySpeciesConcentration,
                                                   logPrimarySpeciesConcentration0,
                                                   logPrimarySpeciesConcentration,
                                                   logSecondarySpeciesConcentration );
}


template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_CONSTRAINT,
          typename ARRAY_1D_SECONDARY,
          typename ARRAY_1D_AGGREGATE >
HPCREACT_HOST_DEVICE inline
bool
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE,
                      ACTIVITY_MODEL >::enforceEquilibrium_PrimaryConcentrations( REAL_TYPE const & temperature,
                                                                                  PARAMS_DATA const & params,
                                                                                  typename ACTIVITY_MODEL::Params const & activityParams,
                                                                                  ARRAY_1D_CONSTRAINT const & constraintType,
                                                                                  ARRAY_1D_TO_CONST const & constraintValue,
                                                                                  ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration0,
                                                                                  ARRAY_1D & logPrimarySpeciesConcentration,
                                                                                  ARRAY_1D_SECONDARY & logSecondarySpeciesConcentration,
                                                                                  ARRAY_1D_AGGREGATE & aggregatePrimarySpeciesConcentration )
{
  HPCREACT_UNUSED_VAR( temperature );
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  double residual[numPrimarySpecies] = { 0.0 };
  double dLogCp[numPrimarySpecies] = { 0.0 };
  CArrayWrapper< double, numPrimarySpecies, numPrimarySpecies > jacobian;

  for( int i=0; i<numPrimarySpecies; ++i )
  {
    logPrimarySpeciesConcentration[i] = logPrimarySpeciesConcentration0[i];
  }

#if HPCREACT_IDEAL_PRESOLVE
  // Generate the initial guess for the full self-consistent solve from the ideal solution of the
  // same system: the one with every activity coefficient and the water activity fixed at 1, which
  // this loop reaches from any starting point.
  //
  // An arbitrary guess, such as the target aggregatePrimaryConcentrations concentrations, can give unrealistic ionic
  // strength, water activity and activity coefficients. The secondary species iteration then stops
  // converging and this loop reaches NaN within a few iterations.
  using IdealActivityModel = Identity< RealType, IndexType, typename ACTIVITY_MODEL::IonicStrengthType >;
  if constexpr( !std::is_same< ACTIVITY_MODEL, IdealActivityModel >::value )
  {
    using IdealEquilibriumReactions = EquilibriumReactions< REAL_TYPE, INT_TYPE, INDEX_TYPE, IdealActivityModel >;

    // Identity reads nothing but the ionic strength parameters it shares with the real model.
    typename IdealActivityModel::Params const idealActivityParams
    {
      static_cast< typename ACTIVITY_MODEL::IonicStrengthType::Params const & >( activityParams )
    };

    IdealEquilibriumReactions::enforceEquilibrium_PrimaryConcentrations( temperature,
                                                                         params,
                                                                         idealActivityParams,
                                                                         constraintType,
                                                                         constraintValue,
                                                                         logPrimarySpeciesConcentration0,
                                                                         logPrimarySpeciesConcentration,
                                                                         logSecondarySpeciesConcentration,
                                                                         aggregatePrimarySpeciesConcentration );
  }
#endif

  // TODO: find an appropriate scaler for the residual.
  constexpr REAL_TYPE residualNormTolerance = 1.0e-8;

  REAL_TYPE residualNorm = 0.0;
  bool isConverged = false;
  bool speciationConverged = true;

  for( int k=0; k<150; ++k )
  {
    speciationConverged &=
      computeResidualAndJacobianPrimaryConcentrations( temperature,
                                                       params,
                                                       activityParams,
                                                       constraintType,
                                                       constraintValue,
                                                       logPrimarySpeciesConcentration,
                                                       residual,
                                                       jacobian,
                                                       logSecondarySpeciesConcentration,
                                                       aggregatePrimarySpeciesConcentration );

    residualNorm = 0.0;
    for( int i = 0; i < numPrimarySpecies; ++i )
    {
      residualNorm += residual[i] * residual[i];
    }
    residualNorm = sqrt( residualNorm );

#if HPCREACT_SOLVER_DIAGNOSTICS
    printf( "iter, residualNorm = %2d, %16.10g \n", k, residualNorm );
#endif

    if( residualNorm < residualNormTolerance )
    {
#if HPCREACT_SOLVER_DIAGNOSTICS
      printf( " converged\n" );
#endif
      isConverged = true;
      break;
    }

    solveNxN_pivoted< double, numPrimarySpecies >( jacobian.data, residual, dLogCp );


    for( IndexType i=0; i<numPrimarySpecies; ++i )
    {
      logPrimarySpeciesConcentration[i] += dLogCp[i];
    }

  }

  return isConverged && speciationConverged;
}

} // namespace reactionsSystems
} // namespace hpcReact
