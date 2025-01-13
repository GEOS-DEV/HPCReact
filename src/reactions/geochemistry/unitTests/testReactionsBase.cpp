
#include "../ReactionsBase_impl.hpp"
#include "MultiVector.hpp"
#include "common/CArrayWrapper.hpp"
#include "../GeochemicalReactionsParametersPredefined.hpp"

#include <gtest/gtest.h>

using namespace hpcReact;

TEST( testReactionsBase, testParamsInitialization )
{
  static_assert( chemicalReactionsParams.numPrimarySpecies == 7, "Number of primary species is not 7" );
  static_assert( chemicalReactionsParams.numSecondarySpecies == 11, "Number of secondary species is not 11" );
  static_assert( chemicalReactionsParams.m_ionSizePrimary[0] == 9.00, "Ion size for primary species 0 is not 9.00" );
  static_assert( chemicalReactionsParams.m_ionSizePrimary[1] == 4.00, "Ion size for primary species 1 is not 4.00" );
  static_assert( chemicalReactionsParams.m_ionSizePrimary[2] == 6.00, "Ion size for primary species 2 is not 6.00" );
  static_assert( chemicalReactionsParams.m_ionSizePrimary[3] == 4.00, "Ion size for primary species 3 is not 4.00" );
  static_assert( chemicalReactionsParams.m_ionSizePrimary[4] == 3.00, "Ion size for primary species 4 is not 3.00" );
  static_assert( chemicalReactionsParams.m_ionSizePrimary[5] == 8.00, "Ion size for primary species 5 is not 8.00" );
  static_assert( chemicalReactionsParams.m_ionSizePrimary[6] == 4.00, "Ion size for primary species 6 is not 4.00" );
  static_assert( chemicalReactionsParams.m_ionSizeSec[0] == 3.50, "Ion size for secondary species 0 is not 3.50" );
  static_assert( chemicalReactionsParams.m_ionSizeSec[1] == 3.00, "Ion size for secondary species 1 is not 3.00" );
  static_assert( chemicalReactionsParams.m_ionSizeSec[2] == 4.50, "Ion size for secondary species 2 is not 4.50" );
  static_assert( chemicalReactionsParams.m_ionSizeSec[3] == 3.00, "Ion size for secondary species 3 is not 3.00" );
  static_assert( chemicalReactionsParams.m_ionSizeSec[4] == 4.00, "Ion size for secondary species 4 is not 4.00" );
  static_assert( chemicalReactionsParams.m_ionSizeSec[5] == 3.00, "Ion size for secondary species 5 is not 3.00" );
  static_assert( chemicalReactionsParams.m_ionSizeSec[6] == 3.00, "Ion size for secondary species 6 is not 3.00" );
  static_assert( chemicalReactionsParams.m_ionSizeSec[7] == 4.00, "Ion size for secondary species 7 is not 4.00" );
  static_assert( chemicalReactionsParams.m_ionSizeSec[8] == 3.00, "Ion size for secondary species 8 is not 3.00" );
  static_assert( chemicalReactionsParams.m_ionSizeSec[9] == 3.00, "Ion size for secondary species 9 is not 3.00" );
  static_assert( chemicalReactionsParams.m_ionSizeSec[10] == 4.00, "Ion size for secondary species 10 is not 4.00" );
  static_assert( chemicalReactionsParams.m_chargePrimary[0] == 1, "Charge for primary species 0 is not 1" );
  static_assert( chemicalReactionsParams.m_chargePrimary[1] == -1, "Charge for primary species 1 is not -1" );
  static_assert( chemicalReactionsParams.m_chargePrimary[2] == 2, "Charge for primary species 2 is not 2" );
  static_assert( chemicalReactionsParams.m_chargePrimary[3] == -2, "Charge for primary species 3 is not -2" );
  static_assert( chemicalReactionsParams.m_chargePrimary[4] == -1, "Charge for primary species 4 is not -1" );
  static_assert( chemicalReactionsParams.m_chargePrimary[5] == 2, "Charge for primary species 5 is not 2" );
  static_assert( chemicalReactionsParams.m_chargePrimary[6] == 1, "Charge for primary species 6 is not 1" );
  static_assert( chemicalReactionsParams.m_chargeSec[0] == -1, "Charge for secondary species 0 is not -1" );
  static_assert( chemicalReactionsParams.m_chargeSec[1] == 0, "Charge for secondary species 1 is not 0" );
  static_assert( chemicalReactionsParams.m_chargeSec[2] == -2, "Charge for secondary species 2 is not -2" );
  static_assert( chemicalReactionsParams.m_chargeSec[3] == 0, "Charge for secondary species 3 is not 0" );
  static_assert( chemicalReactionsParams.m_chargeSec[4] == 1, "Charge for secondary species 4 is not 1" );
  static_assert( chemicalReactionsParams.m_chargeSec[5] == 0, "Charge for secondary species 5 is not 0" );
  static_assert( chemicalReactionsParams.m_chargeSec[6] == 0, "Charge for secondary species 6 is not 0" );
  static_assert( chemicalReactionsParams.m_chargeSec[7] == 1, "Charge for secondary species 7 is not 1" );
  static_assert( chemicalReactionsParams.m_chargeSec[8] == 0, "Charge for secondary species 8 is not 0" );
  static_assert( chemicalReactionsParams.m_chargeSec[9] == 0, "Charge for secondary species 9 is not 0" );
  static_assert( chemicalReactionsParams.m_chargeSec[10] == -1, "Charge for secondary species 10 is not -1" );
  static_assert( ( chemicalReactionsParams.m_DebyeHuckelA - 0.5465 ) * ( chemicalReactionsParams.m_DebyeHuckelA - 0.5465 ) < 1.0e-15, "Debye Huckel A is not 0.5465" );
  static_assert( ( chemicalReactionsParams.m_DebyeHuckelB - 0.3346 ) * ( chemicalReactionsParams.m_DebyeHuckelB - 0.3346 ) < 1.0e-15, "Debye Huckel B is not 0.3346" );
  static_assert( ( chemicalReactionsParams.m_WATEQBDot - 0.0438 ) * ( chemicalReactionsParams.m_WATEQBDot - 0.0438 ) < 1.0e-15, "WATEQBDot is not 0.0438" );
}



template< typename REACTIONS_BASE_TYPE, typename PARAMS_TYPE >
auto testComputIonicStrengthHelper( PARAMS_TYPE const & params,
                                    typename REACTIONS_BASE_TYPE::RealConstDataArrayView1d const & primarySpeciesConcentration,
                                    typename REACTIONS_BASE_TYPE::RealConstDataArrayView1d const & secondarySpeciesConcentration )
{
  using ReactionBaseType = REACTIONS_BASE_TYPE;

  typename ReactionBaseType::RealType ionicStrength = 0.0;

  ReactionBaseType::computeIonicStrength( params,
                                          primarySpeciesConcentration,
                                          secondarySpeciesConcentration,
                                          ionicStrength );

  return ionicStrength;

}

TEST( testReactionsBase, testComputeIonicStrength )
{

  using ReactionsType = ReactionsBase< double,
                                       double *,
                                       double const *,
                                       int,
                                       int >;


  double const primarySpeciesConcentration[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
  double const secondarySpeciesConcentration[] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 };

  typename ReactionsType::RealType const ionicStrength = testComputIonicStrengthHelper< ReactionsType >( chemicalReactionsParams,
                                                                                                         primarySpeciesConcentration,
                                                                                                         secondarySpeciesConcentration );

  EXPECT_NEAR( ionicStrength, 52.0, 1.0e-12 );
}


TEST( testReactionsBase, testComputeIonicStrengthVector )
{

  using ReactionsType = ReactionsBase< double,
                                       std::vector< double >,
                                       std::vector< double const >,
                                       int,
                                       int >;

  std::vector< double const > const primarySpeciesConcentration{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0 };
  std::vector< double const > const secondarySpeciesConcentration{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 };
  typename ReactionsType::RealType const ionicStrength = testComputIonicStrengthHelper< ReactionsType >( chemicalReactionsParams,
                                                                                                         primarySpeciesConcentration,
                                                                                                         secondarySpeciesConcentration );

  EXPECT_NEAR( ionicStrength, 52.0, 1.0e-12 );

}

template< typename REACTIONS_BASE_TYPE, typename PARAMS_TYPE >
void testComputeLog10ActCoefBDotModelHelper( PARAMS_TYPE const & params,
                                             typename REACTIONS_BASE_TYPE::RealType const temperature,
                                             typename REACTIONS_BASE_TYPE::RealType const ionicStrength,
                                             typename REACTIONS_BASE_TYPE::RealDataArrayView1d & log10PrimaryActCoeff,
                                             typename REACTIONS_BASE_TYPE::RealDataArrayView1d & dLog10PrimaryActCoeff_dIonicStrength,
                                             typename REACTIONS_BASE_TYPE::RealDataArrayView1d & log10SecActCoeff,
                                             typename REACTIONS_BASE_TYPE::RealDataArrayView1d & dLog10SecActCoeff_dIonicStrength )
{
  using ReactionBaseType = REACTIONS_BASE_TYPE;

  ReactionBaseType::computeLog10ActCoefBDotModel( temperature,
                                                  ionicStrength,
                                                  params,
                                                  log10PrimaryActCoeff,
                                                  dLog10PrimaryActCoeff_dIonicStrength,
                                                  log10SecActCoeff,
                                                  dLog10SecActCoeff_dIonicStrength );
}

TEST( testReactionsBase, testComputeLog10ActCoefBDotModel )
{
  using ReactionsType = ReactionsBase< double,
                                       double * const,
                                       double const * const,
                                       int,
                                       int >;

  double const temperature = 298.15;
  double const ionicStrength = 0.5;

  double log10PrimaryActCoeff[7] = {0};
  double dLog10PrimaryActCoeff_dIonicStrength[7] = {0};

  double log10SecActCoeff[11] = {0};
  double dLog10SecActCoeff_dIonicStrength[11] = {0};


  testComputeLog10ActCoefBDotModelHelper< ReactionsType >( chemicalReactionsParams,
                                                           temperature,
                                                           ionicStrength,
                                                           log10PrimaryActCoeff,
                                                           dLog10PrimaryActCoeff_dIonicStrength,
                                                           log10SecActCoeff,
                                                           dLog10SecActCoeff_dIonicStrength );

}

TEST( testReactionsBase, testComputeLog10ActCoefBDotModel_vector )
{
  using ReactionsType = ReactionsBase< double,
                                       std::vector< double >,
                                       std::vector< double const >,
                                       int,
                                       int >;

  double const temperature = 298.15;
  double const ionicStrength = 0.5;

  std::vector< double > log10PrimaryActCoeff( 7 );
  std::vector< double > dLog10PrimaryActCoeff_dIonicStrength( 7 );

  std::vector< double > log10SecActCoeff( 11 );
  std::vector< double > dLog10SecActCoeff_dIonicStrength( 11 );


  testComputeLog10ActCoefBDotModelHelper< ReactionsType >( chemicalReactionsParams,
                                                           temperature,
                                                           ionicStrength,
                                                           log10PrimaryActCoeff,
                                                           dLog10PrimaryActCoeff_dIonicStrength,
                                                           log10SecActCoeff,
                                                           dLog10SecActCoeff_dIonicStrength );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
