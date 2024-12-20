
#include "../ReactionsBase_impl.hpp"
#include <gtest/gtest.h>

using namespace hpcReact;

TEST( testReactionsBase, testParamsInitialization )
{
  constexpr
  ReactionsBase< double, 
                 int, 
                 double, 
                 int, 
                 int >::ParamsData< 7, 11 > params { { 9.00, 4.00, 6.00, 4.00, 3.00, 8.00, 4.00 },
                                                     { 3.50, 3.00, 4.50, 3.00, 4.00, 3.00, 3.00, 4.00, 3.00, 3.00, 4.00 },
                                                     { 1, -1, 2, -2, -1, 2, 1 },
                                                     { -1, 0, -2, 0, 1, 0, 0, 1, 0, 0, -1}, 
                                                     0.5465, 
                                                     0.3346, 
                                                     0.0438 };
  static_assert( params.numPrimarySpecies() == 7, "Number of primary species is not 7" );
  static_assert( params.numSecondarySpecies() == 11, "Number of secondary species is not 11" );
  static_assert( params.m_ionSizePrimary[0] == 9.00, "Ion size for primary species 0 is not 9.00" );
  static_assert( params.m_ionSizePrimary[1] == 4.00, "Ion size for primary species 1 is not 4.00" );
  static_assert( params.m_ionSizePrimary[2] == 6.00, "Ion size for primary species 2 is not 6.00" );
  static_assert( params.m_ionSizePrimary[3] == 4.00, "Ion size for primary species 3 is not 4.00" );
  static_assert( params.m_ionSizePrimary[4] == 3.00, "Ion size for primary species 4 is not 3.00" );
  static_assert( params.m_ionSizePrimary[5] == 8.00, "Ion size for primary species 5 is not 8.00" );
  static_assert( params.m_ionSizePrimary[6] == 4.00, "Ion size for primary species 6 is not 4.00" );
  static_assert( params.m_ionSizeSec[0] == 3.50, "Ion size for secondary species 0 is not 3.50" );
  static_assert( params.m_ionSizeSec[1] == 3.00, "Ion size for secondary species 1 is not 3.00" );
  static_assert( params.m_ionSizeSec[2] == 4.50, "Ion size for secondary species 2 is not 4.50" );
  static_assert( params.m_ionSizeSec[3] == 3.00, "Ion size for secondary species 3 is not 3.00" );
  static_assert( params.m_ionSizeSec[4] == 4.00, "Ion size for secondary species 4 is not 4.00" );
  static_assert( params.m_ionSizeSec[5] == 3.00, "Ion size for secondary species 5 is not 3.00" );
  static_assert( params.m_ionSizeSec[6] == 3.00, "Ion size for secondary species 6 is not 3.00" );
  static_assert( params.m_ionSizeSec[7] == 4.00, "Ion size for secondary species 7 is not 4.00" );
  static_assert( params.m_ionSizeSec[8] == 3.00, "Ion size for secondary species 8 is not 3.00" );
  static_assert( params.m_ionSizeSec[9] == 3.00, "Ion size for secondary species 9 is not 3.00" );
  static_assert( params.m_ionSizeSec[10] == 4.00, "Ion size for secondary species 10 is not 4.00" );
  static_assert( params.m_chargePrimary[0] == 1, "Charge for primary species 0 is not 1" );
  static_assert( params.m_chargePrimary[1] == -1, "Charge for primary species 1 is not -1" );
  static_assert( params.m_chargePrimary[2] == 2, "Charge for primary species 2 is not 2" );
  static_assert( params.m_chargePrimary[3] == -2, "Charge for primary species 3 is not -2" );
  static_assert( params.m_chargePrimary[4] == -1, "Charge for primary species 4 is not -1" );
  static_assert( params.m_chargePrimary[5] == 2, "Charge for primary species 5 is not 2" );
  static_assert( params.m_chargePrimary[6] == 1, "Charge for primary species 6 is not 1" );
  static_assert( params.m_chargeSec[0] == -1, "Charge for secondary species 0 is not -1" );
  static_assert( params.m_chargeSec[1] == 0, "Charge for secondary species 1 is not 0" );
  static_assert( params.m_chargeSec[2] == -2, "Charge for secondary species 2 is not -2" );
  static_assert( params.m_chargeSec[3] == 0, "Charge for secondary species 3 is not 0" );
  static_assert( params.m_chargeSec[4] == 1, "Charge for secondary species 4 is not 1" );
  static_assert( params.m_chargeSec[5] == 0, "Charge for secondary species 5 is not 0" );
  static_assert( params.m_chargeSec[6] == 0, "Charge for secondary species 6 is not 0" );
  static_assert( params.m_chargeSec[7] == 1, "Charge for secondary species 7 is not 1" );
  static_assert( params.m_chargeSec[8] == 0, "Charge for secondary species 8 is not 0" );
  static_assert( params.m_chargeSec[9] == 0, "Charge for secondary species 9 is not 0" );
  static_assert( params.m_chargeSec[10] == -1, "Charge for secondary species 10 is not -1" );
  static_assert( ( params.m_DebyeHuckelA - 0.5465 ) * ( params.m_DebyeHuckelA - 0.5465 ) < 1.0e-15, "Debye Huckel A is not 0.5465" );
  static_assert( ( params.m_DebyeHuckelB - 0.3346 ) * ( params.m_DebyeHuckelB - 0.3346 ) < 1.0e-15, "Debye Huckel B is not 0.3346" );
  static_assert( ( params.m_WATEQBDot - 0.0438 ) * ( params.m_WATEQBDot - 0.0438 ) < 1.0e-15, "WATEQBDot is not 0.0438" );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
