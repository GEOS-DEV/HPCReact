
#include "reactions/reactionsSystems/EquilibriumReactions.hpp"
#include "common/macros.hpp"

#include <gtest/gtest.h>


namespace hpcReact
{

namespace unitTest_utilities
{

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
                                     REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies()],
                                     REAL_TYPE const (&expectedResidual)[PARAMS_DATA::numReactions()],
                                     REAL_TYPE const (&expectedJacobian)[PARAMS_DATA::numReactions()][PARAMS_DATA::numReactions()] )
{
  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< REAL_TYPE,
                                                                           int,
                                                                           int >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies();
  constexpr int numReactions = PARAMS_DATA::numReactions();

  double const temperature = 298.15;
  double speciesConcentration[numSpecies];

  for( int i = 0; i < numSpecies; ++i )
  {
    speciesConcentration[i] = initialSpeciesConcentration[i];
  }

  double residual[numReactions] = { 0.0 };
  double xi[numReactions] = { 0.0 };

  CArrayWrapper< double, numReactions, numReactions > jacobian;


  EquilibriumReactionsType::computeResidualAndJacobianReactionExtents( temperature,
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
template< typename REAL_TYPE,
          int RESIDUAL_FORM,
          typename PARAMS_DATA >
void testEnforceEquilibrium( PARAMS_DATA const & params,
                             REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numSpecies()],
                             REAL_TYPE const (&expectedSpeciesConcentrations)[PARAMS_DATA::numSpecies()] )
{
  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< REAL_TYPE,
                                                                           int,
                                                                           int >;

  constexpr int numSpecies = PARAMS_DATA::numSpecies();

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

}

} // namespace hpcReact
