#include "../SpeciesUtilities.hpp"
#include "../ParametersPredefined.hpp"
#include "common/printers.hpp"

#include <gtest/gtest.h>

using namespace hpcReact;
using namespace hpcReact::bulkGeneric;



TEST( testUtilities, test_calculateLogSecondarySpeciesConcentration )
{
  constexpr int numReactions = carbonateSystem.numReactions;
  constexpr int numSpecies = carbonateSystem.numSpecies;
  constexpr int numPrimarySpecies = numSpecies - numReactions;
  constexpr int numSecondarySpecies = numReactions;

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

  double logSecondarySpeciesConcentrations[numSecondarySpecies] = {0};

  calculateLogSecondarySpeciesConcentration< double,
                                             int,
                                             int >( carbonateSystem,
                                                    logPrimarySpeciesSolution,
                                                    logSecondarySpeciesConcentrations );

  double expectedSecondarySpeciesConcentrations[numSecondarySpecies] =
  {
    2.3278416955869804e-11,
    0.3745973700632361,
    3.956656978189425e-11,
    0.0009629355924566718,
    0.00006739226982791149,
    1.065032288527949e-9,
    0.005298329882666744,
    0.005844517547638335,
    0.012773193926706526,
    0.006618125707964999,
    0.01769217213462983
  };

  for( int j=0; j<numSecondarySpecies; ++j )
  {
    EXPECT_NEAR( exp( logSecondarySpeciesConcentrations[j] ),
                 expectedSecondarySpeciesConcentrations[j],
                 1.0e-8 );
  }


  double dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[numSecondarySpecies][numPrimarySpecies] = {{0}};

  calculateLogSecondarySpeciesConcentrationWrtLogC< double,
                                                    int,
                                                    int >( carbonateSystem,
                                                           logPrimarySpeciesSolution,
                                                           logSecondarySpeciesConcentrations,
                                                           dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations );



  double expected_dLogSdLogC[numSecondarySpecies][numPrimarySpecies] =
  {
    { -1, 0, 0, 0, 0, 0, 0 },
    {  1, 1, 0, 0, 0, 0, 0 },
    { -1, 1, 0, 0, 0, 0, 0 },
    {  1, 1, 0, 0, 0, 0, 0 },
    {  0, 1, 1, 0, 0, 0, 0 },
    { -1, 1, 1, 0, 0, 0, 0 },
    {  0, 0, 1, 1, 0, 0, 0 },
    {  0, 0, 1, 0, 1, 0, 0 },
    {  0, 0, 1, 0, 2, 0, 0 },
    {  0, 0, 0, 1, 0, 1, 0 },
    {  0, 0, 0, 1, 0, 0, 1 }
  };

  for( int i=0; i<numSecondarySpecies; ++i )
  {
    for( int j=0; j<numPrimarySpecies; ++j )
    {
      EXPECT_NEAR( dLogSecondarySpeciesConcentrations_dLogPrimarySpeciesConcentrations[i][j],
                   expected_dLogSdLogC[i][j],
                   1.0e-8 );
    }
  }
}


TEST( testUtilities, testcalculateAggregatePrimaryConcentrationsWrtLogC )
{
  constexpr int numReactions = carbonateSystem.numReactions;
  constexpr int numSpecies = carbonateSystem.numSpecies;
  constexpr int numPrimarySpecies = numSpecies - numReactions;
//  constexpr int numSecondarySpecies = numReactions;

  double primarySpeciesSolution[numPrimarySpecies] =
  {
    0.00043969547214915125,
    0.00037230096984514874,
    0.014716565308128551,
    0.0024913722747387217,
    1.8586090945989489,
    0.009881874292035079,
    1.0723078278653704
  };

  for( int i=0; i<numPrimarySpecies; ++i )
  {
    primarySpeciesSolution[i] = log( primarySpeciesSolution[i] );
  }


  double aggregatePrimarySpeciesConcentration[numPrimarySpecies] = {0};


  CArrayWrapper< double, numPrimarySpecies, numPrimarySpecies > dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations = {{{0.0}}};

  calculateAggregatePrimaryConcentrationsWrtLogC< double, int, int >( carbonateSystem,
                                                                      primarySpeciesSolution,
                                                                      aggregatePrimarySpeciesConcentration,
                                                                      dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations );

  double const expectedAggregatePrimarySpeciesConcentration[numPrimarySpecies] =
  {
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89, // Cl-
    1.65e-2, // Mg+2
    1.09 // Na+1
  };

  for( int i=0; i<numPrimarySpecies; ++i )
  {
    EXPECT_NEAR( aggregatePrimarySpeciesConcentration[i],
                 expectedAggregatePrimarySpeciesConcentration[i],
                 1.0e-8 );
  }

  double expected_dAggregatePrimarySpeciesConcentration_dLogPrimarySpeciesConcentrations[numPrimarySpecies][numPrimarySpecies] =
  {
    { 0.3760000022557546, 0.3755603045511293, -1.0650322885265095e-9, 0, 0, 0, 0 },
    { 0.3755603045511293, 0.37600000000000006, 0.00006739333486015734, 0, 0, 0, 0 },
    { -1.0650322885265095e-9, 0.00006739333486015734, 0.03870000000000001, 0.005298329882666651, 0.031390905401051036, 0, 0 },
    { 0, 0, 0.005298329882666651, 0.032100000000000004, 0, 0.006618125707964927, 0.01769217213462971 },
    { 0, 0, 0.031390905401051036, 0, 1.9155463878534127, 0, 0 },
    { 0, 0, 0, 0.006618125707964927, 0, 0.016500000000000008, 0 },
    { 0, 0, 0, 0.01769217213462971, 0, 0, 1.09 }
  };


  for( int i=0; i<numPrimarySpecies; ++i )
  {
    for( int j=0; j<numPrimarySpecies; ++j )
    {
      EXPECT_NEAR( dAggregatePrimarySpeciesConcentrationsDerivatives_dLogPrimarySpeciesConcentrations( i, j ),
                   expected_dAggregatePrimarySpeciesConcentration_dLogPrimarySpeciesConcentrations[i][j],
                   1.0e-8 );
    }
  }
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
