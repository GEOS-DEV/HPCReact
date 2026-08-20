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

#include "../MassActions.hpp"
#include "reactions/geochemistry/GeochemicalSystems.hpp"
#include "common/pmpl.hpp"
#include "common/printers.hpp"
#include <math.h>


#include <gtest/gtest.h>

using namespace hpcReact;
using namespace hpcReact::massActions;
using namespace hpcReact::geochemistry;


template< int numPrimarySpecies, int numSecondarySpecies >
struct CalculateLogSecondarySpeciesConcentrationNoActivityUpdateData
{
  double logSecondarySpeciesConcentrations[numSecondarySpecies] = {0};
  double dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[numSecondarySpecies][numPrimarySpecies] = {{0}};
  double const logPrimarySpeciesSolution[numPrimarySpecies] =
  {
    log( 0.00043969547214915125 ),
    log( 0.00037230096984514874 ),
    log( 0.014716565308128551 ),
    log( 0.0024913722747387217 ),
    log( 1.8586090945989489 ),
    log( 0.009881874292035079 ),
    log( 1.0723078278653704 )
  };
};


void test_calculateLogSecondarySpeciesConcentrationNoActivityUpdate_helper()
{
  static constexpr int numPrimarySpecies = carbonateSystemAllEquilibrium.numPrimarySpecies();
  static constexpr int numSecondarySpecies = carbonateSystemAllEquilibrium.numSecondarySpecies();

  CalculateLogSecondarySpeciesConcentrationNoActivityUpdateData< numPrimarySpecies, numSecondarySpecies > data;

  // What the Identity model would return: gamma = 1, so the secondary activities and
  // concentrations coincide.
  double const identityLogSecondaryActivityCoefficients[numSecondarySpecies] = { 0.0 };

  pmpl::genericKernelWrapper( 1, &data, [identityLogSecondaryActivityCoefficients] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    calculateLogSecondarySpeciesConcentrationNoActivityUpdate< double,
                                                               int,
                                                               int >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                                      dataCopy->logPrimarySpeciesSolution,
                                                                      identityLogSecondaryActivityCoefficients,
                                                                      dataCopy->logSecondarySpeciesConcentrations );



    calculateLogSecondarySpeciesConcentrationWrtLogCNoActivityUpdate< double,
                                                                      int,
                                                                      int >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                                             dataCopy->logPrimarySpeciesSolution,
                                                                             identityLogSecondaryActivityCoefficients,
                                                                             dataCopy->logSecondarySpeciesConcentrations,
                                                                             dataCopy->dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations );
  } );



  double expectedSecondarySpeciesConcentrations[numSecondarySpecies] =
  {
    2.3278416955869804e-11,
    0.37035984325260107,
    3.956656978189425e-11,
    9.1316525616762169e-05,
    0.007654372189572687,
    0.13676171061473555,
    0.012773193926706526,
    0.0041586871002752598,
    0.013225336595688551,
    0.00024148987937519677
  };

  for( int j=0; j<numSecondarySpecies; ++j )
  {
    EXPECT_NEAR( exp( data.logSecondarySpeciesConcentrations[j] ),
                 expectedSecondarySpeciesConcentrations[j],
                 1.0e-8 );
  }

  double expected_dLogSdLogC[numSecondarySpecies][numPrimarySpecies] =
  {
    { -1, 0, 0, 0, 0, 0, 0 },
    {  1, 1, 0, 0, 0, 0, 0 },
    { -1, 1, 0, 0, 0, 0, 0 },
    {  0, 1, 1, 0, 0, 0, 0 },
    {  0, 0, 1, 1, 0, 0, 0 },
    {  0, 0, 1, 0, 1, 0, 0 },
    {  0, 0, 1, 0, 2, 0, 0 },
    {  0, 0, 0, 1, 0, 1, 0 },
    {  0, 0, 0, 1, 0, 0, 1 },
    { -1, 1, 1, 0, 0, 0, 0 }

  };

  for( int i=0; i<numSecondarySpecies; ++i )
  {
    for( int j=0; j<numPrimarySpecies; ++j )
    {
      EXPECT_NEAR( data.dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[i][j],
                   expected_dLogSdLogC[i][j],
                   1.0e-8 );
    }
  }
}

TEST( testMassActions, test_calculateLogSecondarySpeciesConcentrationNoActivityUpdate )
{
  test_calculateLogSecondarySpeciesConcentrationNoActivityUpdate_helper();
}


template< int numPrimarySpecies, int numSecondarySpecies, int numSpecies >
struct CalculateLogSecondarySpeciesConcentrationData
{
  double logSecondarySpeciesConcentrations[numSecondarySpecies] = {0};
  double logActivityCoefficients[numSpecies] = {0};
  double logActivities[numSpecies] = {0};
  double dLogActivities_dLogSpeciesConcentrations[numSpecies][numSpecies] = {{0}};
  double dLogActivityCoefficients_dLogSpeciesConcentrations[numSpecies][numSpecies] = {{0}};
  bool converged = false;

  double const logPrimarySpeciesSolution[numPrimarySpecies] =
  {
    log( 0.00043969547214915125 ),
    log( 0.00037230096984514874 ),
    log( 0.014716565308128551 ),
    log( 0.0024913722747387217 ),
    log( 1.8586090945989489 ),
    log( 0.009881874292035079 ),
    log( 1.0723078278653704 )
  };
};

/**
 * @brief Verify the coupled secondary concentration and activity coefficient solve is
 *   self-consistent.
 * @param massActionTolerance tolerance on the equilibrium residual.
 *
 * Ionic strength couples the secondary concentrations and the activity coefficients, so a Newton
 * iteration updates both together. The checks below are on the returned state itself, which is why
 * one helper serves any activity model:
 *   - mass action holds on the activities, to massActionTolerance
 *   - a = C*gamma for the secondary species, to 1e-12
 *   - a = C*gamma for the primary species, to 1e-12
 */
template< typename ACTIVITY_MODEL, typename EQ_PARAMS, typename ACTIVITY_PARAMS >
void test_calculateLogSecondarySpeciesConcentration_helper( EQ_PARAMS const eqParams,
                                                            ACTIVITY_PARAMS const activityParams,
                                                            double const massActionTolerance )
{
  static constexpr int numPrimarySpecies   = EQ_PARAMS::numPrimarySpecies();
  static constexpr int numSecondarySpecies = EQ_PARAMS::numSecondarySpecies();
  static constexpr int numSpecies          = EQ_PARAMS::numSpecies();

  CalculateLogSecondarySpeciesConcentrationData< numPrimarySpecies, numSecondarySpecies, numSpecies > data;

  pmpl::genericKernelWrapper( 1, &data, [eqParams, activityParams] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    dataCopy->converged =
      calculateLogSecondarySpeciesConcentration< double,
                                                 int,
                                                 int,
                                                 ACTIVITY_MODEL,
                                                 true >( eqParams,
                                                         activityParams,
                                                         dataCopy->logPrimarySpeciesSolution,
                                                         dataCopy->logSecondarySpeciesConcentrations,
                                                         dataCopy->logActivityCoefficients,
                                                         dataCopy->logActivities,
                                                         dataCopy->dLogActivities_dLogSpeciesConcentrations,
                                                         dataCopy->dLogActivityCoefficients_dLogSpeciesConcentrations );
  } );

  EXPECT_TRUE( data.converged );

  // The equilibrium constraint must hold at the returned state:
  //   log(a_j) + log(K_j) - sum_k nu_jk log(a_k) = 0
  for( int j = 0; j < numSecondarySpecies; ++j )
  {
    double residual = data.logActivities[j] + log( eqParams.equilibriumConstant( j ) );
    for( int k = 0; k < numPrimarySpecies; ++k )
    {
      residual -= eqParams.stoichiometricMatrix( j, k + numSecondarySpecies ) *
                  data.logActivities[k + numSecondarySpecies];
    }
    EXPECT_NEAR( residual, 0.0, massActionTolerance ) << "mass action violated, secondary species " << j;
  }

  // Concentrations, activity coefficients and activities must be mutually consistent.
  for( int j = 0; j < numSecondarySpecies; ++j )
  {
    EXPECT_NEAR( data.logActivities[j],
                 data.logSecondarySpeciesConcentrations[j] + data.logActivityCoefficients[j],
                 1.0e-12 ) << "a != c*gamma, secondary species " << j;
  }

  // The primary concentrations were held fixed, so their activities follow from them directly.
  for( int k = 0; k < numPrimarySpecies; ++k )
  {
    EXPECT_NEAR( data.logActivities[k + numSecondarySpecies],
                 data.logPrimarySpeciesSolution[k] + data.logActivityCoefficients[k + numSecondarySpecies],
                 1.0e-12 ) << "a != c*gamma, primary species " << k;
  }
}

// Two systems, each with the Identity and B-dot activity models: carbonateSystemAllEquilibrium
// (all reactions at equilibrium) and carbonateSystem (mixed equilibrium and kinetic).

TEST( testMassActions, test_calculateLogSecondarySpeciesConcentration_allEquilibrium_identity )
{
  test_calculateLogSecondarySpeciesConcentration_helper< carbonateIdentityActivityType >(
    carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
    carbonateIdentityActivityParams,
    1.0e-12 );
}

TEST( testMassActions, test_calculateLogSecondarySpeciesConcentration_allEquilibrium_bdot )
{
  test_calculateLogSecondarySpeciesConcentration_helper< carbonateActivityType >(
    carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
    carbonateActivityParamsEQ36,
    1.0e-9 );
}

TEST( testMassActions, test_calculateLogSecondarySpeciesConcentration_mixedSystem_identity )
{
  test_calculateLogSecondarySpeciesConcentration_helper< carbonateNosolidIdentityActivityType >(
    carbonateSystem.equilibriumReactionsParameters(),
    carbonateNosolidIdentityActivityParams,
    1.0e-12 );
}

TEST( testMassActions, test_calculateLogSecondarySpeciesConcentration_mixedSystem_bdot )
{
  test_calculateLogSecondarySpeciesConcentration_helper< carbonateNosolidActivityType >(
    carbonateSystem.equilibriumReactionsParameters(),
    carbonateNosolidActivityParams,
    1.0e-9 );
}

/**
 * @brief With the Identity model the solve must reproduce the reference solution.
 */
void test_calculateLogSecondarySpeciesConcentration_identityActivityModel_helper()
{
  static constexpr int numPrimarySpecies = carbonateSystemAllEquilibrium.numPrimarySpecies();
  static constexpr int numSecondarySpecies = carbonateSystemAllEquilibrium.numSecondarySpecies();
  static constexpr int numSpecies = carbonateSystemAllEquilibrium.numSpecies();

  CalculateLogSecondarySpeciesConcentrationData< numPrimarySpecies, numSecondarySpecies, numSpecies > data;

  constexpr auto activityParams = carbonateIdentityActivityParams;

  pmpl::genericKernelWrapper( 1, &data, [activityParams] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    dataCopy->converged =
      calculateLogSecondarySpeciesConcentration< double,
                                                 int,
                                                 int,
                                                 carbonateIdentityActivityType,
                                                 true >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                         activityParams,
                                                         dataCopy->logPrimarySpeciesSolution,
                                                         dataCopy->logSecondarySpeciesConcentrations,
                                                         dataCopy->logActivityCoefficients,
                                                         dataCopy->logActivities,
                                                         dataCopy->dLogActivities_dLogSpeciesConcentrations,
                                                         dataCopy->dLogActivityCoefficients_dLogSpeciesConcentrations );
  } );

  // Reference solution, as in test_calculateLogSecondarySpeciesConcentrationNoActivityUpdate above:
  // gamma = 1 makes the ideal seed exact, so the Newton has nothing to correct.
  double const expectedSecondarySpeciesConcentrations[numSecondarySpecies] =
  {
    2.3278416955869804e-11,
    0.37035984325260107,
    3.956656978189425e-11,
    9.1316525616762169e-05,
    0.007654372189572687,
    0.13676171061473555,
    0.012773193926706526,
    0.0041586871002752598,
    0.013225336595688551,
    0.00024148987937519677
  };

  for( int j = 0; j < numSecondarySpecies; ++j )
  {
    EXPECT_NEAR( exp( data.logSecondarySpeciesConcentrations[j] ),
                 expectedSecondarySpeciesConcentrations[j],
                 1.0e-8 );
    EXPECT_NEAR( data.logActivityCoefficients[j], 0.0, 1.0e-15 );
  }
}

TEST( testMassActions, test_calculateLogSecondarySpeciesConcentration_identityActivityModel )
{
  test_calculateLogSecondarySpeciesConcentration_identityActivityModel_helper();
}


template< int numPrimarySpecies, int numSecondarySpecies, int numSpecies >
struct SpeciationDerivativeData : public CalculateLogSecondarySpeciesConcentrationData< numPrimarySpecies, numSecondarySpecies, numSpecies >
{
  double analytic[numSecondarySpecies][numPrimarySpecies] = {{0}};
  double finiteDifference[numSecondarySpecies][numPrimarySpecies] = {{0}};
};

/**
 * @brief Check d log(C_sec)/d log(C_prim) against a central difference of the solve itself.
 *
 * This derivative feeds the outer Newton that solves for the primary concentrations given the
 * aggregate totals for equilibrium reaction run.
 */
template< typename ACTIVITY_MODEL, typename EQ_PARAMS, typename ACTIVITY_PARAMS >
void test_calculateLogSecondarySpeciesConcentrationWrtLogC_helper( EQ_PARAMS const eqParams,
                                                                   ACTIVITY_PARAMS const activityParams,
                                                                   double const tolerance )
{
  static constexpr int numPrimarySpecies   = EQ_PARAMS::numPrimarySpecies();
  static constexpr int numSecondarySpecies = EQ_PARAMS::numSecondarySpecies();
  static constexpr int numSpecies          = EQ_PARAMS::numSpecies();

  SpeciationDerivativeData< numPrimarySpecies, numSecondarySpecies, numSpecies > data;

  pmpl::genericKernelWrapper( 1, &data, [eqParams, activityParams] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    dataCopy->converged =
      calculateLogSecondarySpeciesConcentrationWrtLogC< double,
                                                        int,
                                                        int,
                                                        ACTIVITY_MODEL,
                                                        true >( eqParams,
                                                                activityParams,
                                                                dataCopy->logPrimarySpeciesSolution,
                                                                dataCopy->logSecondarySpeciesConcentrations,
                                                                dataCopy->logActivityCoefficients,
                                                                dataCopy->logActivities,
                                                                dataCopy->dLogActivities_dLogSpeciesConcentrations,
                                                                dataCopy->dLogActivityCoefficients_dLogSpeciesConcentrations,
                                                                dataCopy->analytic );

    double constexpr perturbation = 1.0e-6;

    for( int n = 0; n < numPrimarySpecies; ++n )
    {
      double logSecondaryPerturbed[2][numSecondarySpecies] = {{0}};

      for( int side = 0; side < 2; ++side )
      {
        double logPrimaryPerturbed[numPrimarySpecies] = {0};
        for( int k = 0; k < numPrimarySpecies; ++k )
        {
          logPrimaryPerturbed[k] = dataCopy->logPrimarySpeciesSolution[k];
        }
        logPrimaryPerturbed[n] += ( side == 0 ? perturbation : -perturbation );

        // Scratch: only the concentrations are wanted here.
        double logActivityCoefficients[numSpecies] = {0};
        double logActivities[numSpecies] = {0};
        double dLogActivities_dLogSpeciesConcentrations[numSpecies][numSpecies] = {{0}};
        double dLogActivityCoefficients_dLogSpeciesConcentrations[numSpecies][numSpecies] = {{0}};

        calculateLogSecondarySpeciesConcentration< double,
                                                   int,
                                                   int,
                                                   ACTIVITY_MODEL,
                                                   true >( eqParams,
                                                           activityParams,
                                                           logPrimaryPerturbed,
                                                           logSecondaryPerturbed[side],
                                                           logActivityCoefficients,
                                                           logActivities,
                                                           dLogActivities_dLogSpeciesConcentrations,
                                                           dLogActivityCoefficients_dLogSpeciesConcentrations );
      }

      for( int j = 0; j < numSecondarySpecies; ++j )
      {
        dataCopy->finiteDifference[j][n] =
          ( logSecondaryPerturbed[0][j] - logSecondaryPerturbed[1][j] ) / ( 2.0 * perturbation );
      }
    }
  } );

  EXPECT_TRUE( data.converged );

  for( int j = 0; j < numSecondarySpecies; ++j )
  {
    for( int n = 0; n < numPrimarySpecies; ++n )
    {
      EXPECT_NEAR( data.analytic[j][n], data.finiteDifference[j][n], tolerance )
        << "d log(C_sec[" << j << "])/d log(C_prim[" << n << "])";
    }
  }
}

TEST( testMassActions, test_calculateLogSecondarySpeciesConcentrationWrtLogC_allEquilibrium_identity )
{
  test_calculateLogSecondarySpeciesConcentrationWrtLogC_helper< carbonateIdentityActivityType >(
    carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
    carbonateIdentityActivityParams,
    1.0e-7 );
}

TEST( testMassActions, test_calculateLogSecondarySpeciesConcentrationWrtLogC_allEquilibrium_bdot )
{
  test_calculateLogSecondarySpeciesConcentrationWrtLogC_helper< carbonateActivityType >(
    carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
    carbonateActivityParamsEQ36,
    1.0e-7 );
}

TEST( testMassActions, test_calculateLogSecondarySpeciesConcentrationWrtLogC_mixedSystem_bdot )
{
  test_calculateLogSecondarySpeciesConcentrationWrtLogC_helper< carbonateNosolidActivityType >(
    carbonateSystem.equilibriumReactionsParameters(),
    carbonateNosolidActivityParams,
    1.0e-7 );
}


template< int numPrimarySpecies, int numSecondarySpecies, int numSpecies >
struct CalculateAggregatePrimaryConcentrationsWrtLogCHelperData
  : public CalculateLogSecondarySpeciesConcentrationData< numPrimarySpecies, numSecondarySpecies, numSpecies >
{
  double dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[numSecondarySpecies][numPrimarySpecies] = {{0}};
  double aggregatePrimarySpeciesConcentration[numPrimarySpecies] = {0};
  CArrayWrapper< double, numPrimarySpecies, numPrimarySpecies > dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations;

};

/**
 * @brief Check the aggregate mole balance against the ideal-solution reference values.
 *
 * The aggregate is a pure mole balance on whatever secondary state it is handed, so the reference
 * values below are reached through the Identity activity model rather than through a separate
 * ideal-solution entry point.
 */
void testcalculateAggregatePrimaryConcentrationsWrtLogCHelper()
{
  static constexpr int numPrimarySpecies   = carbonateSystemAllEquilibrium.numPrimarySpecies();
  static constexpr int numSecondarySpecies = carbonateSystemAllEquilibrium.numSecondarySpecies();
  static constexpr int numSpecies          = carbonateSystemAllEquilibrium.numSpecies();

  CalculateAggregatePrimaryConcentrationsWrtLogCHelperData< numPrimarySpecies, numSecondarySpecies, numSpecies > data;

  constexpr auto activityParams = carbonateIdentityActivityParams;

  pmpl::genericKernelWrapper( 1, &data, [activityParams] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    dataCopy->converged =
      calculateLogSecondarySpeciesConcentrationWrtLogC< double,
                                                        int,
                                                        int,
                                                        carbonateIdentityActivityType,
                                                        true >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                                activityParams,
                                                                dataCopy->logPrimarySpeciesSolution,
                                                                dataCopy->logSecondarySpeciesConcentrations,
                                                                dataCopy->logActivityCoefficients,
                                                                dataCopy->logActivities,
                                                                dataCopy->dLogActivities_dLogSpeciesConcentrations,
                                                                dataCopy->dLogActivityCoefficients_dLogSpeciesConcentrations,
                                                                dataCopy->dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations );

    calculateAggregatePrimaryConcentrationsWrtLogC< double, int, int >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                                        dataCopy->logPrimarySpeciesSolution,
                                                                        dataCopy->logSecondarySpeciesConcentrations,
                                                                        dataCopy->dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations,
                                                                        dataCopy->aggregatePrimarySpeciesConcentration,
                                                                        dataCopy->dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
  } );

  EXPECT_TRUE( data.converged );

  double const expectedAggregatePrimarySpeciesConcentration[numPrimarySpecies] =
  {
    0.37055804878406567, // H+
    0.37106495066575157, // HCO3-
    0.17223864844413511, // Ca+2
    0.02752976816027522, // SO4-2
    2.0209171930670973, // Cl-
    0.014040561392310339, // Mg+2
    1.0855331644610589 // Na+1
  };

  for( int i=0; i<numPrimarySpecies; ++i )
  {
    EXPECT_NEAR( data.aggregatePrimarySpeciesConcentration[i],
                 expectedAggregatePrimarySpeciesConcentration[i],
                 1.0e-8 );
  }

  double expected_dAggregatePrimarySpeciesConcentration_dLogPrimarySpeciesConcentrations[numPrimarySpecies][numPrimarySpecies] =
  {
    { 0.37104102866543476, 0.3701183533349125, -0.00024148987937519677, 0, 0, 0, 0 },
    { 0.3701183533349125, 0.37106495066575157, 0.00033280640499195897, 0, 0, 0, 0 },
    { -0.00024148987937519677, 0.00033280640499195897, 0.17223864844413511, 0.007654372189572687, 0.16230809846814831, 0, 0 },
    { 0, 0, 0.007654372189572687, 0.02752976816027522, 0, 0.0041586871002752598, 0.013225336595688551 },
    { 0, 0, 0.16230809846814831, 0, 2.0464635809205101, 0, 0 },
    { 0, 0, 0, 0.0041586871002752598, 0, 0.014040561392310339, 0 },
    { 0, 0, 0, 0.013225336595688551, 0, 0, 1.0855331644610589 }
  };


  for( int i=0; i<numPrimarySpecies; ++i )
  {
    for( int j=0; j<numPrimarySpecies; ++j )
    {
      EXPECT_NEAR( data.dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, j ),
                   expected_dAggregatePrimarySpeciesConcentration_dLogPrimarySpeciesConcentrations[i][j],
                   1.0e-8 );
    }
  }
}

TEST( testMassActions, testcalculateAggregatePrimaryConcentrationsWrtLogC )
{
  testcalculateAggregatePrimaryConcentrationsWrtLogCHelper();
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
