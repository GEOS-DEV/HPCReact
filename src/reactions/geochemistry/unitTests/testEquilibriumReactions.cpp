
#include "../EquilibriumReactions_impl.hpp"

#include <gtest/gtest.h>

using namespace hpcReact;

TEST( testEquilibriumReactions, setInitialGuess )
{}


TEST( testEquilibriumReactions, assembleEquilibriumReactionSystem )
{}


TEST( testEquilibriumReactions, computeSecondarySpeciesConcAndDerivative )
{}


TEST( testEquilibriumReactions, computeTotalConcAndDerivative )
{}

TEST( testEquilibriumReactions, updatePrimarySpeciesConcentrations )
{}

TEST( testEquilibriumReactions, updateConcentrations )
{
  
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
