
#include "../EquilibriumReactions.hpp"
#include "../ParametersPredefined.hpp"
#include "common/macros.hpp"

#include <gtest/gtest.h>


using namespace hpcReact;
using namespace hpcReact::bulkGeneric;


template< typename REAL_TYPE >
REAL_TYPE tolerance( REAL_TYPE const a, REAL_TYPE const b )
{
  return std::numeric_limits< double >::epsilon() * std::max( fabs( a ), fabs( b ) ) * 10;
}

//******************************************************************************
template< typename REAL_TYPE,
          int RESIDUAL_FORM,
          typename PARAMS_DATA >
void computeResidualAndJacobianTest( PARAMS_DATA const & params,
                                     REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies],
                                     REAL_TYPE const (&expectedResidual)[PARAMS_DATA::numReactions],
                                     REAL_TYPE const (&expectedJacobian)[PARAMS_DATA::numReactions][PARAMS_DATA::numReactions] )
{
  using EquilibriumReactionsType = EquilibriumReactions< REAL_TYPE,
                                                         int,
                                                         int >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies;
  constexpr int numReactions = PARAMS_DATA::numReactions;

  double const temperature = 298.15;
  double speciesConcentration[numSpecies];

  for( int i = 0; i < numSpecies; ++i )
  {
    speciesConcentration[i] = initialSpeciesConcentration[i];
  }

  double residual[numReactions] = { 0.0 };
  double xi[numReactions] = { 0.0 };

  CArrayWrapper< double, numReactions, numReactions > jacobian;


  EquilibriumReactionsType::computeResidualAndJacobian( temperature,
                                                        params,
                                                        speciesConcentration,
                                                        xi,
                                                        residual,
                                                        jacobian );

//  printf( "R = { %8.4g, %8.4g }\n", residual[0], residual[1] );
  for( int r=0; r<numReactions; ++r )
  {
    EXPECT_NEAR( residual[r], expectedResidual[r], 1.0e-8 );
  }

  // HPCREACT_UNUSED_VAR( expectedJacobian );
  // printf( "J = { { %8.4g, %8.4g },\n    { %8.4g %8.4g } }\n\n",
  //         jacobian( 0, 0 ), jacobian( 0, 1 ),
  //         jacobian( 1, 0 ), jacobian( 1, 1 ) );
  for( int a = 0; a < numReactions; ++a )
  {
    for( int b = 0; b < numReactions; ++b )
    {
      EXPECT_NEAR( jacobian( a, b ), expectedJacobian[a][b], tolerance( jacobian( a, b ), expectedJacobian[a][b] ) );
    }
  }
}


//******************************************************************************
TEST( testEquilibriumReactions, computeResidualAndJacobianTest )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };


  {
    std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
    double const expectedResiduals[] = { -37.534508668465, -72.989575795250 };
    double const expectedJacobian[2][2] =
    { { 1.0e16, -2.0 },
      { -2.0, 4.0e16 } };
    computeResidualAndJacobianTest< double, 2 >( simpleTestRateParams,
                                                 initialSpeciesConcentration,
                                                 expectedResiduals,
                                                 expectedJacobian );
  }

}


//******************************************************************************
template< typename REAL_TYPE,
          int RESIDUAL_FORM,
          typename PARAMS_DATA >
void testEnforceEquilibrium( PARAMS_DATA const & params,
                             REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies],
                             REAL_TYPE const (&expectedSpeciesConcentrations)[PARAMS_DATA::numSpecies] )
{
  using EquilibriumReactionsType = EquilibriumReactions< REAL_TYPE,
                                                         int,
                                                         int >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies;

  double const temperature = 298.15;
  double speciesConcentration0[numSpecies];
  double speciesConcentration[numSpecies];


  for( int i = 0; i < numSpecies; ++i )
  {
    speciesConcentration0[i] = initialSpeciesConcentration[i];
  }

  EquilibriumReactionsType::enforceEquilibrium_Extents( temperature,
                                                params,
                                                speciesConcentration0,
                                                speciesConcentration );

  for( int r=0; r<numSpecies; ++r )
  {
//    printf( "c[%d] = %22.14e\n", r, speciesConcentration[r] );
    EXPECT_NEAR( speciesConcentration[r], expectedSpeciesConcentrations[r], 1.0e-8 * expectedSpeciesConcentrations[r] );
  }

}


//******************************************************************************
TEST( testEquilibriumReactions, testEnforceEquilibrium )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedSpeciesConcentrations[5] = { 3.92138294e-01, 3.03930853e-01, 5.05945481e-01, 7.02014628e-01, 5.95970745e-01 };


  std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
  testEnforceEquilibrium< double, 2 >( simpleTestRateParams.equilibriumReactionsParameters(),
                                       initialSpeciesConcentration,
                                       expectedSpeciesConcentrations );

}


//******************************************************************************
TEST( testEquilibriumReactions, testCarbonateSystem )
{
  double const initialSpeciesConcentration[18] =
  {
    1.0e-16, // OH-
    1.0e-16, // CO2
    1.0e-16, // CO3-2
    1.0e-16, // H2CO3
    1.0e-16, // CaHCO3+
    1.0e-16, // CaCO3
    1.0e-16, // CaSO4
    1.0e-16, // CaCl+
    1.0e-16, // CaCl2
    1.0e-16, // MgSO4
    1.0e-16, // NaSO4-
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89, // Cl-
    1.65e-2, // Mg+2
    1.09 // Na+1
  };

  double const expectedSpeciesConcentrations[18] =
  { 2.327841695586879e-11, // OH-
    3.745973700632716e-01, // CO2
    3.956656978189456e-11, // CO3-2
    9.629355924567627e-04, // H2CO3
    6.739226982791492e-05, // CaHCO3+
    1.065032288527957e-09, // CaCO3
    5.298329882666738e-03, // CaSO4
    5.844517547638333e-03, // CaCl+
    1.277319392670652e-02, // CaCl2
    6.618125707964991e-03, // MgSO4
    1.769217213462983e-02, // NaSO4-
    4.396954721488358e-04, // H+
    3.723009698453808e-04, // HCO3-
    1.471656530812871e-02, // Ca+2
    2.491372274738741e-03, // SO4-2
    1.858609094598949e+00, // Cl-
    9.881874292035110e-03, // Mg+2
    1.072307827865370e+00 // Na+1
  };

  std::cout<<" RESIDUAL_FORM 0:"<<std::endl;
  testEnforceEquilibrium< double, 0 >( carbonateSystem.equilibriumReactionsParameters(),
                                       initialSpeciesConcentration,
                                       expectedSpeciesConcentrations );

  // std::cout<<" RESIDUAL_FORM 1:"<<std::endl;
  // testEnforceEquilibrium< double, 1 >( carbonateSystem.equilibriumReactionsParameters(),
  //                                      initialSpeciesConcentration,
  //                                      expectedSpeciesConcentrations );

  std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
  testEnforceEquilibrium< double, 2 >( carbonateSystem.equilibriumReactionsParameters(),
                                       initialSpeciesConcentration,
                                       expectedSpeciesConcentrations );

}


TEST( testEquilibriumReactions, testCarbonateSystem2 )
{
  using EquilibriumReactionsType = EquilibriumReactions< double,
                                                         int,
                                                         int >;

  constexpr int numSpecies = carbonateSystem.numSpecies;
  constexpr int numReactions = carbonateSystem.numReactions;
  constexpr int numPrimarySpecies = numSpecies - numReactions;
//  constexpr int numSecondarySpecies = numReactions;

  double const initialPrimarySpeciesConcentration[numPrimarySpecies] =
  {
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89,    // Cl-
    1.65e-2, // Mg+2
    1.09     // Na+1
  };



  double const logInitialPrimarySpeciesConcentration[numPrimarySpecies] = 
  {
    log( initialPrimarySpeciesConcentration[0] ),
    log( initialPrimarySpeciesConcentration[1] ),
    log( initialPrimarySpeciesConcentration[2] ),
    log( initialPrimarySpeciesConcentration[3] ),
    log( initialPrimarySpeciesConcentration[4] ),
    log( initialPrimarySpeciesConcentration[5] ),
    log( initialPrimarySpeciesConcentration[6] )
  };

  double logPrimarySpeciesConcentration[numPrimarySpecies] = 
  {
    logInitialPrimarySpeciesConcentration[0],
    logInitialPrimarySpeciesConcentration[1],
    logInitialPrimarySpeciesConcentration[2],
    logInitialPrimarySpeciesConcentration[3],
    logInitialPrimarySpeciesConcentration[4],
    logInitialPrimarySpeciesConcentration[5],
    logInitialPrimarySpeciesConcentration[6]
  };

  EquilibriumReactionsType::enforceEquilibrium_Aggregate( 0,
                                                          carbonateSystem,
                                                          logInitialPrimarySpeciesConcentration,
                                                          logPrimarySpeciesConcentration );

  double const expectedPrimarySpeciesConcentrations[numPrimarySpecies] =
  { 
    4.396954721488358e-04, // H+
    3.723009698453808e-04, // HCO3-
    1.471656530812871e-02, // Ca+2
    2.491372274738741e-03, // SO4-2
    1.858609094598949e+00, // Cl-
    9.881874292035110e-03, // Mg+2
    1.072307827865370e+00 // Na+1
  };

  for( int r=0; r<numPrimarySpecies; ++r )
  {
    EXPECT_NEAR( exp( logPrimarySpeciesConcentration[r] ), expectedPrimarySpeciesConcentrations[r], 1.0e-8 * expectedPrimarySpeciesConcentrations[r] );
  }


}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
