
#include "../EquilibriumReactions.hpp"
#include "../ParametersPredefined.hpp"
#include "common/DirectSystemSolve.hpp"
#include "common/macros.hpp"

#include <gtest/gtest.h>


using namespace hpcReact;
using namespace hpcReact::bulkGeneric;


//******************************************************************************
template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void computeResidualAndJacobianTest( PARAMS_DATA const & params,
                                     REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies],
                                     REAL_TYPE const (&expectedResidual)[PARAMS_DATA::numReactions],
                                     REAL_TYPE const (&expectedJacobian)[PARAMS_DATA::numReactions][PARAMS_DATA::numReactions] )
{
  using EquilibriumReactionsType = EquilibriumReactions< REAL_TYPE,
                                                         int,
                                                         int,
                                                         false,
                                                         false >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies;
  constexpr int numReactions = PARAMS_DATA::numReactions;

  double const temperature = 298.15;
  double speciesConcentration[numSpecies];

  if constexpr( LOGE_CONCENTRATION )
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      speciesConcentration[i] = log( initialSpeciesConcentration[i] );
    }
  }
  else
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      speciesConcentration[i] = initialSpeciesConcentration[i];
    }
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

  for( int r=0; r<numReactions; ++r )
  {
    EXPECT_NEAR( residual[r], expectedResidual[r], 1.0e-8 );
  }

  for( int a = 0; a < numReactions; ++a )
  {
    for( int b = 0; b < numReactions; ++b )
    {
      EXPECT_NEAR( jacobian( a, b ), expectedJacobian[a][b], 1.0e-8 );
    }
  }
}


//******************************************************************************
TEST( testEquilibriumReactions, computeResidualAndJacobianTest )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedResiduals[] = { 1.0, 0.5 };
  double const expectedJacobian[2][2] = { { -4.5, 0 },
    { 1.0, -1.5 } };
  computeResidualAndJacobianTest< double, false >( simpleTestRateParams,
                                                   initialSpeciesConcentration,
                                                   expectedResiduals,
                                                   expectedJacobian );
}


//******************************************************************************
template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void testEnforceEquilibrium( PARAMS_DATA const & params,
                             REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies],
                             REAL_TYPE const (&expectedSpeciesConcentrations)[PARAMS_DATA::numSpecies] )
{
  using EquilibriumReactionsType = EquilibriumReactions< REAL_TYPE,
                                                         int,
                                                         int,
                                                         false,
                                                         false >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies;

  double const temperature = 298.15;
  double speciesConcentration[numSpecies];

  if constexpr( LOGE_CONCENTRATION )
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      speciesConcentration[i] = log( initialSpeciesConcentration[i] );
    }
  }
  else
  {
    for( int i = 0; i < numSpecies; ++i )
    {
      speciesConcentration[i] = initialSpeciesConcentration[i];
    }
  }

  EquilibriumReactionsType::enforceEquilibrium( temperature,
                                                params,
                                                initialSpeciesConcentration,
                                                speciesConcentration );

  for( int r=0; r<numSpecies; ++r )
  {
    EXPECT_NEAR( speciesConcentration[r], expectedSpeciesConcentrations[r], 1.0e-8 );
  }

}


//******************************************************************************
TEST( testEquilibriumReactions, testEnforceEquilibrium )
{
  double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
  double const expectedSpeciesConcentrations[5] = { 3.92138294e-01, 3.03930853e-01, 5.05945481e-01, 7.02014628e-01, 5.95970745e-01 };
  testEnforceEquilibrium< double, false >( simpleTestRateParams,
                                           initialSpeciesConcentration,
                                           expectedSpeciesConcentrations );
}



int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
