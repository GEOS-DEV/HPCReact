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


#include <gtest/gtest.h>

using namespace hpcReact;
using namespace hpcReact::massActions;
using namespace hpcReact::geochemistry;


template< int numPrimarySpecies, int numSecondarySpecies >
struct CalculateLogSecondarySpeciesConcentrationData
{
  double logSecondarySpeciesConcentrations[numSecondarySpecies] = {0};
  double dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[numSecondarySpecies][numPrimarySpecies] = {{0}};
  double const logPrimarySpeciesSolution[numPrimarySpecies] =
  {
    logmath::log( 0.00043969547214915125 ),
    logmath::log( 0.00037230096984514874 ),
    logmath::log( 0.014716565308128551 ),
    logmath::log( 0.0024913722747387217 ),
    logmath::log( 1.8586090945989489 ),
    logmath::log( 0.009881874292035079 ),
    logmath::log( 1.0723078278653704 )
  };
};


void test_calculateLogSecondarySpeciesConcentration_helper()
{
  static constexpr int numPrimarySpecies = carbonateSystemAllEquilibrium.numPrimarySpecies();
  static constexpr int numSecondarySpecies = carbonateSystemAllEquilibrium.numSecondarySpecies();

  CalculateLogSecondarySpeciesConcentrationData< numPrimarySpecies, numSecondarySpecies > data;

  pmpl::genericKernelWrapper( 1, &data, [] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    calculateLogSecondarySpeciesConcentration< double,
                                               int,
                                               int >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                      dataCopy->logPrimarySpeciesSolution,
                                                      dataCopy->logSecondarySpeciesConcentrations );



    calculateLogSecondarySpeciesConcentrationWrtLogC< double,
                                                      int,
                                                      int >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                             dataCopy->logPrimarySpeciesSolution,
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
    EXPECT_NEAR( logmath::exp( data.logSecondarySpeciesConcentrations[j] ),
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

TEST( testMassActions, test_calculateLogSecondarySpeciesConcentration )
{
  test_calculateLogSecondarySpeciesConcentration_helper();
}


template< int numPrimarySpecies >
struct CalculateAggregatePrimaryConcentrationsWrtLogCHelperData
{
  double const primarySpeciesSolution[numPrimarySpecies] =
  {
    logmath::log( 0.00043969547214915125 ),
    logmath::log( 0.00037230096984514874 ),
    logmath::log( 0.014716565308128551 ),
    logmath::log( 0.0024913722747387217 ),
    logmath::log( 1.8586090945989489 ),
    logmath::log( 0.009881874292035079 ),
    logmath::log( 1.0723078278653704 )
  };

  double aggregatePrimarySpeciesConcentration[numPrimarySpecies] = {0};
  CArrayWrapper< double, numPrimarySpecies, numPrimarySpecies > dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations;

};

void testcalculateAggregatePrimaryConcentrationsWrtLogCHelper()
{
  static constexpr int numPrimarySpecies = carbonateSystemAllEquilibrium.numPrimarySpecies();

  CalculateAggregatePrimaryConcentrationsWrtLogCHelperData< numPrimarySpecies > data;

  pmpl::genericKernelWrapper( 1, &data, [] HPCREACT_DEVICE ( auto * const dataCopy )
  {
    calculateAggregatePrimaryConcentrationsWrtLogC< double, int, int >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                                        dataCopy->primarySpeciesSolution,
                                                                        dataCopy->aggregatePrimarySpeciesConcentration,
                                                                        dataCopy->dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );
  } );

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
                   logmath::dWrtLogConst<double>() * expected_dAggregatePrimarySpeciesConcentration_dLogPrimarySpeciesConcentrations[i][j],
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
