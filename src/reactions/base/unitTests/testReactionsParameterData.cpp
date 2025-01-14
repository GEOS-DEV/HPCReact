
#include "../ReactionsBase_impl.hpp"
#include "MultiVector.hpp"
#include "common/CArrayWrapper.hpp"
#include "../ParametersPredefined.hpp"

#include <gtest/gtest.h>

using namespace hpcReact;

template< typename T >
constexpr bool nearEqual( T a, T b, T tol = 1e-7 )
{
  return ( a - b ) * ( a - b ) < tol * tol;
}

TEST( testReactionsBase, testParamsInitialization )
{
  static_assert( chemicalReactionsParams.numPrimarySpecies == 7, "Number of primary species is not 7" );
  static_assert( chemicalReactionsParams.numSecondarySpecies == 11, "Number of secondary species is not 11" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizePrimary[0], 9.00 ), "Ion size for primary species 0 is not 9.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizePrimary[1], 4.00 ), "Ion size for primary species 1 is not 4.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizePrimary[2], 6.00 ), "Ion size for primary species 2 is not 6.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizePrimary[3], 4.00 ), "Ion size for primary species 3 is not 4.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizePrimary[4], 3.00 ), "Ion size for primary species 4 is not 3.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizePrimary[5], 8.00 ), "Ion size for primary species 5 is not 8.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizePrimary[6], 4.00 ), "Ion size for primary species 6 is not 4.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizeSec[0] , 3.50 ), "Ion size for secondary species 0 is not 3.50" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizeSec[1] , 3.00 ), "Ion size for secondary species 1 is not 3.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizeSec[2] , 4.50 ), "Ion size for secondary species 2 is not 4.50" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizeSec[3] , 3.00 ), "Ion size for secondary species 3 is not 3.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizeSec[4] , 4.00 ), "Ion size for secondary species 4 is not 4.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizeSec[5] , 3.00 ), "Ion size for secondary species 5 is not 3.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizeSec[6] , 3.00 ), "Ion size for secondary species 6 is not 3.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizeSec[7] , 4.00 ), "Ion size for secondary species 7 is not 4.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizeSec[8] , 3.00 ), "Ion size for secondary species 8 is not 3.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizeSec[9] , 3.00 ), "Ion size for secondary species 9 is not 3.00" );
  static_assert( nearEqual( chemicalReactionsParams.m_ionSizeSec[10], 4.00 ), "Ion size for secondary species 10 is not 4.00" );
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
  static_assert( nearEqual( chemicalReactionsParams.m_DebyeHuckelA, 0.5465 ), "Debye Huckel A is not 0.5465" );
  static_assert( nearEqual( chemicalReactionsParams.m_DebyeHuckelB, 0.3346 ), "Debye Huckel B is not 0.3346" );
  static_assert( nearEqual( chemicalReactionsParams.m_WATEQBDot, 0.0438 ), "WATEQBDot is not 0.0438" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[0][0], -1.0 ), "Stoichiometry matrix for species 0 and 0 is not -1" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[0][1],  0.0 ), "Stoichiometry matrix for species 0 and 1 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[0][2],  0.0 ), "Stoichiometry matrix for species 0 and 2 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[0][3],  0.0 ), "Stoichiometry matrix for species 0 and 3 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[0][4],  0.0 ), "Stoichiometry matrix for species 0 and 4 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[0][5],  0.0 ), "Stoichiometry matrix for species 0 and 5 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[0][6],  0.0 ), "Stoichiometry matrix for species 0 and 6 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[1][0],  1.0 ), "Stoichiometry matrix for species 1 and 0 is not 1" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[1][1],  1.0 ), "Stoichiometry matrix for species 1 and 1 is not 1" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[1][2],  0.0 ), "Stoichiometry matrix for species 1 and 2 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[1][3],  0.0 ), "Stoichiometry matrix for species 1 and 3 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[1][4],  0.0 ), "Stoichiometry matrix for species 1 and 4 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[1][5],  0.0 ), "Stoichiometry matrix for species 1 and 5 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqStoichMatrix[1][6],  0.0 ), "Stoichiometry matrix for species 1 and 6 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqLog10EqConst[0], 13.99 ),  "Log10 equilibrium constant for species 0 is not 13.99" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqLog10EqConst[1], -6.36 ) , "Log10 equilibrium constant for species 1 is not -6.36" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqLog10EqConst[2], 10.33 ),  "Log10 equilibrium constant for species 2 is not 10.33" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqLog10EqConst[3], -3.77 ) , "Log10 equilibrium constant for species 3 is not -3.77" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqLog10EqConst[4], -1.09 ) , "Log10 equilibrium constant for species 4 is not -1.09" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqLog10EqConst[5],  7.07 ) , "Log10 equilibrium constant for species 5 is not 7.07" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqLog10EqConst[6], - 2.16 ) , "Log10 equilibrium constant for species 6 is not -2.16" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqLog10EqConst[7], + 0.67 ) , "Log10 equilibrium constant for species 7 is not 0.67" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqLog10EqConst[8], + 0.60 ) , "Log10 equilibrium constant for species 8 is not 0.60" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqLog10EqConst[9],  - 2.43 ) , "Log10 equilibrium constant for species 9 is not -2.43" );
  static_assert( nearEqual( chemicalReactionsParams.m_eqLog10EqConst[10], - 0.82 ),  "Log10 equilibrium constant for species 10 is not -0.82" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[0][0],-2.0 ), "Stoichiometry matrix for kinetic reaction 0 and species 0 is not -2" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[0][1], 0.0 ), "Stoichiometry matrix for kinetic reaction 0 and species 1 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[0][2], 1.0 ), "Stoichiometry matrix for kinetic reaction 0 and species 2 is not 1" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[0][3], 0.0 ), "Stoichiometry matrix for kinetic reaction 0 and species 3 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[0][4], 0.0 ), "Stoichiometry matrix for kinetic reaction 0 and species 4 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[0][5], 0.0 ), "Stoichiometry matrix for kinetic reaction 0 and species 5 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[0][6], 0.0 ), "Stoichiometry matrix for kinetic reaction 0 and species 6 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[1][0],-1.0 ), "Stoichiometry matrix for kinetic reaction 1 and species 0 is not -1" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[1][1], 1.0 ), "Stoichiometry matrix for kinetic reaction 1 and species 1 is not 1" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[1][2], 1.0 ), "Stoichiometry matrix for kinetic reaction 1 and species 2 is not 1" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[1][3], 0.0 ), "Stoichiometry matrix for kinetic reaction 1 and species 3 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[1][4], 0.0 ), "Stoichiometry matrix for kinetic reaction 1 and species 4 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[1][5], 0.0 ), "Stoichiometry matrix for kinetic reaction 1 and species 5 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticStoichMatrix[1][6], 0.0 ), "Stoichiometry matrix for kinetic reaction 1 and species 6 is not 0" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticlog10EqConst[0], 20.19 ), "Log10 equilibrium constant for kinetic reaction 0 is not 20.19" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticlog10EqConst[1], 1.32 ), "Log10 equilibrium constant for kinetic reaction 1 is not 1.32" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticReactionRateConstant[0], 9.95e-1 ), "Reaction rate constant for kinetic reaction 0 is not 9.95e-1" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticReactionRateConstant[1], 9.95e-3 ), "Reaction rate constant for kinetic reaction 1 is not 9.95e-3" );
  static_assert( nearEqual( chemicalReactionsParams.m_kineticSpecificSurfaceArea, 1.0 ), "Specific surface area is not 1.0" );


}

TEST( testReactionsBase, testEquilibriumParamsExtraction )
{
  constexpr auto eqParams = chemicalReactionsParams.equilibriumReactions();
  static_assert( eqParams.numPrimarySpecies == 7, "Number of primary species is not 7" );
  static_assert( eqParams.numSecondarySpecies == 11, "Number of secondary species is not 11" );
  static_assert( nearEqual( eqParams.m_ionSizePrimary[0], 9.00 ), "Ion size for primary species 0 is not 9.00" );
  static_assert( nearEqual( eqParams.m_ionSizePrimary[1], 4.00 ), "Ion size for primary species 1 is not 4.00" );
  static_assert( nearEqual( eqParams.m_ionSizePrimary[2], 6.00 ), "Ion size for primary species 2 is not 6.00" );
  static_assert( nearEqual( eqParams.m_ionSizePrimary[3], 4.00 ), "Ion size for primary species 3 is not 4.00" );
  static_assert( nearEqual( eqParams.m_ionSizePrimary[4], 3.00 ), "Ion size for primary species 4 is not 3.00" );
  static_assert( nearEqual( eqParams.m_ionSizePrimary[5], 8.00 ), "Ion size for primary species 5 is not 8.00" );
  static_assert( nearEqual( eqParams.m_ionSizePrimary[6], 4.00 ), "Ion size for primary species 6 is not 4.00" );
  static_assert( nearEqual( eqParams.m_ionSizeSec[0] , 3.50 ), "Ion size for secondary species 0 is not 3.50" );
  static_assert( nearEqual( eqParams.m_ionSizeSec[1] , 3.00 ), "Ion size for secondary species 1 is not 3.00" );
  static_assert( nearEqual( eqParams.m_ionSizeSec[2] , 4.50 ), "Ion size for secondary species 2 is not 4.50" );
  static_assert( nearEqual( eqParams.m_ionSizeSec[3] , 3.00 ), "Ion size for secondary species 3 is not 3.00" );
  static_assert( nearEqual( eqParams.m_ionSizeSec[
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
