
#include "../KineticReactions_impl.hpp"
#include "../ParametersPredefined.hpp"

#include <gtest/gtest.h>


using namespace hpcReact;
using namespace hpcReact::bulkGeneric;

TEST( bulkGeneric, test_computeReactionRates )
{
  using KineticReactionsType = KineticReactions< double, 
                                                 double * const,
                                                 double const * const, 
                                                 int, 
                                                 int >;

  double const temperature = 298.15;
  double const primarySpeciesConcentration[] = { 1.0, 1.0, 1.0, 1.0, 1.0 };
  double reactionRates[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
  
  KineticReactionsType::computeReactionRates( temperature, 
                                              bicarbonateBuffer, 
                                              primarySpeciesConcentration, 
                                              reactionRates );

  EXPECT_NEAR( reactionRates[0], 0.0, 1.0e-8 );
  EXPECT_NEAR( reactionRates[1], 0.0, 1.0e-8 );
  EXPECT_NEAR( reactionRates[2], 0.0, 1.0e-8 );
  EXPECT_NEAR( reactionRates[3], 0.0, 1.0e-8 );
  EXPECT_NEAR( reactionRates[4], 0.0, 1.0e-8 );
  
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
