
#include "../EquilibriumReactions_impl.hpp"

#include <gtest/gtest.h>

using namespace hpcReact;

TEST( testReactionsBase, testParamsInitialization )
{}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
