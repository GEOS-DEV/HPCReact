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

#include "reactions/massActions/MassActions.hpp"
#include "constitutive/activity/Identity.hpp"

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
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D,
          typename ARRAY_1D_SECONDARY >
HPCREACT_HOST_DEVICE
inline
bool
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE,
                      ACTIVITY_MODEL >::computeResidualAndJacobianAggregatePrimaryConcentrations( RealType const & temperature,
                                                                                                  PARAMS_DATA const & params,
                                                                                                  typename ACTIVITY_MODEL::Params const & activityParams,
                                                                                                  ARRAY_1D_TO_CONST const & targetAggregatePrimaryConcentrations,
                                                                                                  ARRAY_1D_TO_CONST2 const & logPrimarySpeciesConcentration,
                                                                                                  ARRAY_1D & residual,
                                                                                                  ARRAY_2D & jacobian,
                                                                                                  ARRAY_1D_SECONDARY & logSecondarySpeciesConcentration )
{
  HPCREACT_UNUSED_VAR( temperature );
  static constexpr int numSpecies = PARAMS_DATA::numSpecies();
  static constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  static constexpr int numSecondarySpeciesStorage = numSecondarySpecies > 0 ? numSecondarySpecies : 1;
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  bool speciationConverged = true;
  RealType aggregatePrimaryConcentrations[numPrimarySpecies] = {0.0};
  RealType dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[numSecondarySpeciesStorage][numPrimarySpecies] = {{0.0}};
  ARRAY_2D dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations = {{{0.0}}};

  if constexpr( numPrimarySpecies > 0 )
  {
    RealType logActivities[numSpecies] = {0.0};
    RealType dLogActivities_dLogSpeciesConcentrations[numSpecies][numSpecies] = {{0.0}};
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
      HPCREACT_UNUSED_VAR( activityParams );
    }
  }

  // Pure mole balance on the solved state: the activity model is already accounted for in the
  // calculation of the two secondary species arrays, so nothing here reconstructs it.
  massActions::calculateAggregatePrimaryConcentrationsWrtLogC< REAL_TYPE, INT_TYPE, INDEX_TYPE >( params,
                                                                                                  logPrimarySpeciesConcentration,
                                                                                                  logSecondarySpeciesConcentration,
                                                                                                  dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                                                                                                  aggregatePrimaryConcentrations,
                                                                                                  dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );


  for( IndexType i=0; i<numPrimarySpecies; ++i )
  {
    residual[i] = -(1.0 - aggregatePrimaryConcentrations[i] / targetAggregatePrimaryConcentrations[i]);
    for( IndexType j=0; j<numPrimarySpecies; ++j )
    {
      jacobian( i, j ) = -dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations[i][j] / targetAggregatePrimaryConcentrations[i];
    }
  }

  return speciationConverged;
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          typename ACTIVITY_MODEL >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST >
HPCREACT_HOST_DEVICE inline
bool
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE,
                      ACTIVITY_MODEL >::enforceEquilibrium_LogAggregate( REAL_TYPE const & temperature,
                                                                         PARAMS_DATA const & params,
                                                                         typename ACTIVITY_MODEL::Params const & activityParams,
                                                                         ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration0,
                                                                         ARRAY_1D & logPrimarySpeciesConcentration )
{
  HPCREACT_UNUSED_VAR( temperature );
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();
  double targetAggregatePrimarySpeciesConcentration[numPrimarySpecies] = { 0.0 };



  for( int i=0; i<numPrimarySpecies; ++i )
  {
    targetAggregatePrimarySpeciesConcentration[i] = exp( logPrimarySpeciesConcentration0[i] );
  }

  return enforceEquilibrium_Aggregate( temperature,
                                       params,
                                       activityParams,
                                       targetAggregatePrimarySpeciesConcentration,
                                       logPrimarySpeciesConcentration0,
                                       logPrimarySpeciesConcentration );
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
                      ACTIVITY_MODEL >::enforceEquilibrium_Aggregate( REAL_TYPE const & temperature,
                                                                      PARAMS_DATA const & params,
                                                                      typename ACTIVITY_MODEL::Params const & activityParams,
                                                                      ARRAY_1D_TO_CONST const & targetAggregatePrimarySpeciesConcentration,
                                                                      ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration0,
                                                                      ARRAY_1D & logPrimarySpeciesConcentration,
                                                                      ARRAY_1D_SECONDARY & logSecondarySpeciesConcentration )
{
  if constexpr( PARAMS_DATA::numSecondarySpecies() <= 0 )
  {
    return true;
  }

  HPCREACT_UNUSED_VAR( temperature );
  static constexpr int numPrimarySpecies = PARAMS_DATA::numPrimarySpecies();

  double residual[numPrimarySpecies] = { 0.0 };
//  double aggregatePrimarySpeciesConcentration[numPrimarySpecies] = { 0.0 };
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
  // An arbitrary guess, such as the target aggregate concentrations, can give unrealistic ionic
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

    IdealEquilibriumReactions::enforceEquilibrium_Aggregate( temperature,
                                                             params,
                                                             idealActivityParams,
                                                             targetAggregatePrimarySpeciesConcentration,
                                                             logPrimarySpeciesConcentration0,
                                                             logPrimarySpeciesConcentration );
  }
#endif

  REAL_TYPE residualNorm = 0.0;
  bool isConverged = false;
  bool speciationConverged = true;

  for( int k=0; k<150; ++k )
  {
    speciationConverged &=
      computeResidualAndJacobianAggregatePrimaryConcentrations( temperature,
                                                                params,
                                                                activityParams,
                                                                targetAggregatePrimarySpeciesConcentration,
                                                                logPrimarySpeciesConcentration,
                                                                residual,
                                                                jacobian,
                                                                logSecondarySpeciesConcentration );

    residualNorm = 0.0;
    for( int i = 0; i < numPrimarySpecies; ++i )
    {
      residualNorm += residual[i] * residual[i];
    }
    residualNorm = sqrt( residualNorm );

    if( residualNorm < 1.0e-12 )
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
