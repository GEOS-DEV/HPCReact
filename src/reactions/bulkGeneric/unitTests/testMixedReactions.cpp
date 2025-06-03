
#include "../MixedEquilibriumKineticReactions.hpp"
#include "../EquilibriumReactions.hpp"
#include "../ParametersPredefined.hpp"
#include "common/macros.hpp"
#include "common/printers.hpp"
#include "common/nonlinearSolvers.hpp"

#include <gtest/gtest.h>


using namespace hpcReact;
using namespace hpcReact::bulkGeneric;

template< typename REAL_TYPE >
REAL_TYPE tolerance( REAL_TYPE const a, REAL_TYPE const b, REAL_TYPE const ndigits )
{
  return std::numeric_limits< double >::epsilon() * std::max( fabs( a ), fabs( b ) ) * pow( 10, ndigits );
}

//******************************************************************************
template< typename REAL_TYPE,
          bool LOGE_CONCENTRATION,
          typename PARAMS_DATA >
void timeStepTest( PARAMS_DATA const & params,
                   REAL_TYPE const dt,
                   int const numSteps,
                   REAL_TYPE const (&initialSpeciesConcentration)[PARAMS_DATA::numPrimarySpecies()],
                   REAL_TYPE const (&expectedSpeciesConcentrations)[PARAMS_DATA::numPrimarySpecies()] )
{
  HPCREACT_UNUSED_VAR( expectedSpeciesConcentrations );

  using MixedReactionsType = MixedEquilibriumKineticReactions< REAL_TYPE,
                                                                  int,
                                                                  int,
                                                               LOGE_CONCENTRATION >;
  using EquilibriumReactionsType = EquilibriumReactions< REAL_TYPE,
                                                         int,
                                                         int >;

  
  // constexpr int numSpecies = PARAMS_DATA::numSpecies();
  constexpr int numPrimarySpecies   = PARAMS_DATA::numPrimarySpecies();
  constexpr int numSecondarySpecies = PARAMS_DATA::numSecondarySpecies();
  constexpr int numKineticReactions = PARAMS_DATA::numKineticReactions();
   
  // define variables
  double const temperature = 298.15;
  REAL_TYPE logPrimarySpeciesConcentration[numPrimarySpecies];
  // must use CArrayWrapper to ensure correct capture in the lambda functions
  CArrayWrapper< REAL_TYPE, numSecondarySpecies > logSecondarySpeciesConcentration;
  CArrayWrapper< REAL_TYPE, numPrimarySpecies > aggregatePrimarySpeciesConcentration;
  CArrayWrapper< REAL_TYPE, numPrimarySpecies > aggregatePrimarySpeciesConcentration_n;
  CArrayWrapper< REAL_TYPE, numPrimarySpecies, numPrimarySpecies > dAggregatePrimarySpeciesConcentrations_dlogPrimarySpeciesConcentration;
  CArrayWrapper< REAL_TYPE, numKineticReactions > reactionRates;
  CArrayWrapper< REAL_TYPE, numKineticReactions, numPrimarySpecies > dReactionRates_dlogPrimarySpeciesConcentration;
  CArrayWrapper< REAL_TYPE, numPrimarySpecies > aggregateSpeciesRates;
  CArrayWrapper< REAL_TYPE, numPrimarySpecies, numPrimarySpecies > dAggregateSpeciesRates_dlogPrimarySpeciesConcentration;
  

  // Initialize species concentrations
  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    logPrimarySpeciesConcentration[i] = ::log( initialSpeciesConcentration[i] );
    aggregatePrimarySpeciesConcentration[i] = initialSpeciesConcentration[i];
  }
  
  EquilibriumReactionsType::enforceEquilibrium_Aggregate( temperature,
                                                          carbonateSystem.equilibriumReactionsParameters(),
                                                          logPrimarySpeciesConcentration,
                                                          logPrimarySpeciesConcentration );
                                                                                                         
  
  /// Time step loop 
  double time = 0.0;
  for( int t = 0; t < numSteps; ++t )
  { 
    printf( "Timestep %d, Time = %.6e\n", t, time );

    for( int i=0; i < numPrimarySpecies; ++i )
    {
      aggregatePrimarySpeciesConcentration_n[i] = aggregatePrimarySpeciesConcentration[i];
    }

    auto computeResidualAndJacobian = [&] HPCREACT_HOST_DEVICE ( REAL_TYPE const (&logPrimarySpeciesConcentration)[numPrimarySpecies], 
                                                                 REAL_TYPE (&r)[numPrimarySpecies],
                                                                 REAL_TYPE (&J)[numPrimarySpecies][numPrimarySpecies] )  
    { 
      MixedReactionsType::updateMixedSystem( temperature,
                                             params,
                                             logPrimarySpeciesConcentration,
                                             logSecondarySpeciesConcentration,
                                             aggregatePrimarySpeciesConcentration,
                                             dAggregatePrimarySpeciesConcentrations_dlogPrimarySpeciesConcentration,
                                             reactionRates,
                                             dReactionRates_dlogPrimarySpeciesConcentration,
                                             aggregateSpeciesRates,
                                             dAggregateSpeciesRates_dlogPrimarySpeciesConcentration );

                                  
      for ( int i = 0; i < numPrimarySpecies; ++i ) 
      {
        r[i] = ( aggregatePrimarySpeciesConcentration[i] - aggregatePrimarySpeciesConcentration_n[i] ) - aggregateSpeciesRates[i] * dt;
        for ( int j = 0; j < numPrimarySpecies; ++j )
        {
          J[i][j] = dAggregatePrimarySpeciesConcentrations_dlogPrimarySpeciesConcentration[i][j] - dAggregateSpeciesRates_dlogPrimarySpeciesConcentration[i][j] * dt;
        }
      } 
    };
    
    nonlinearSolvers::newtonRaphson< numPrimarySpecies >( logPrimarySpeciesConcentration, computeResidualAndJacobian );

    time += dt;
  }
  
  EXPECT_NEAR( reactionRatesDerivatives( r, i ), expectedReactionRatesDerivatives[r][i], std::max( magScale, fabs( expectedReactionRatesDerivatives[r][i] ) ) * 1.0e-8 );
}

TEST( testMixedReactions, testTimeStep_carbonateSystem )
{
  constexpr int numPrimarySpecies = carbonateSystemType::numPrimarySpecies();

  double const initialSpeciesConcentration[numPrimarySpecies] =
  {
    3.76e-3, // CaCO3 
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89,    // Cl-
    1.65e-2, // Mg+2
    1.09     // Na+1
  };

  double const expectedSpeciesConcentrations[numPrimarySpecies] =
  {
    1.0e-3,  // CaCO3
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89,    // Cl-
    1.65e-2, // Mg+2
    1.09     // Na+1
  };

  timeStepTest< double, true >( carbonateSystem,
                                0.2,
                                10,
                                initialSpeciesConcentration,
                                expectedSpeciesConcentrations );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
