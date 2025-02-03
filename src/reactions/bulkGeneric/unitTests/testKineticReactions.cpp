
#include "../KineticReactions_impl.hpp"
#include "../ParametersPredefined.hpp"

#include <gtest/gtest.h>


using namespace hpcReact;
using namespace hpcReact::bulkGeneric;

// TEST( bulkGeneric, test_computeReactionRates )
// {
//   using KineticReactionsType = KineticReactions< double, 
//                                                  double * const,
//                                                  double const * const, 
//                                                  int, 
//                                                  int >;

//   double const temperature = 298.15;
//   double primarySpeciesConcentration[6] = { 0.01, 0.01, 0.01, 0.01, 0.01, 1.0 };
//   double reactionRates[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  
//   KineticReactionsType::computeReactionRates( temperature, 
//                                               bicarbonateBuffer, 
//                                               primarySpeciesConcentration, 
//                                               reactionRates );

//   printf( "reactionRates = {%8.4e, %8.4e, %8.4e, %8.4e, %8.4e, %8.4e }\n", 
//           reactionRates[0],
//           reactionRates[1],
//           reactionRates[2],
//           reactionRates[3],
//           reactionRates[4],
//           reactionRates[5] );
//   // EXPECT_NEAR( reactionRates[0], 0.0, 1.0e-8 );
//   // EXPECT_NEAR( reactionRates[1], 0.0, 1.0e-8 );
//   // EXPECT_NEAR( reactionRates[2], 0.0, 1.0e-8 );
//   // EXPECT_NEAR( reactionRates[3], 0.0, 1.0e-8 );
//   // EXPECT_NEAR( reactionRates[4], 0.0, 1.0e-8 );
  
// }

TEST( bulkGeneric, test_computeReactionRatesIntegral )
{
  using KineticReactionsType = KineticReactions< double, 
                                                 double * const,
                                                 double const * const, 
                                                 int, 
                                                 int >;

  double const temperature = 298.15;
  double primarySpeciesConcentration[] = { 1.0, 0.0, 0.5, 1.0, 0.0 };
  double reactionRates[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

  constexpr double dt = 1.0e-1;
  constexpr int numPrimarySpecies = decltype(simpleTest)::numPrimarySpecies;
//  constexpr int numKineticReactions = decltype(simpleTest)::numKineticReactions;

  //double time = 0.0;
  for( int i = 0; i < 1000; ++i )
  {
    KineticReactionsType::computeReactionRates( temperature, 
                                                simpleTest, 
                                                primarySpeciesConcentration, 
                                                reactionRates );

    for( int j = 0; j < numPrimarySpecies; ++j )
    {
      primarySpeciesConcentration[j] += dt * reactionRates[j];
    }
    //time += dt;
    // if( i % 10 == 0 )
    // {
    //   printf( " {% 12.8e} { % 12.8e, % 12.8e, % 12.8e, % 12.8e, % 12.8e }\n", 
    //           time,
    //           primarySpeciesConcentration[0],
    //           primarySpeciesConcentration[1],
    //           primarySpeciesConcentration[2],
    //           primarySpeciesConcentration[3],
    //           primarySpeciesConcentration[4] );
    // }
  }

  EXPECT_NEAR( primarySpeciesConcentration[0], 3.92138294e-01, 1.0e-8 );
  EXPECT_NEAR( primarySpeciesConcentration[1], 3.03930853e-01, 1.0e-8 );
  EXPECT_NEAR( primarySpeciesConcentration[2], 5.05945481e-01, 1.0e-8 );
  EXPECT_NEAR( primarySpeciesConcentration[3], 7.02014628e-01, 1.0e-8 );
  EXPECT_NEAR( primarySpeciesConcentration[4], 5.95970745e-01, 1.0e-8 );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
