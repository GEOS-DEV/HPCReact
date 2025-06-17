#include "reactions/unitTestUtilities/mixedReactionsTestUtilities.hpp"
#include "../GeochemicalSystems.hpp"


using namespace hpcReact;
using namespace hpcReact::unitTest_utilities;


TEST( testMixedReactions, testTimeStep_carbonateSystem )
{
  using namespace hpcReact::geochemistry;

  constexpr int numPrimarySpecies = carbonateSystemType::numPrimarySpecies();

  double const initialAggregateSpeciesConcentration[numPrimarySpecies] =
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
    3.3318075516669661e-05,  // CaCO3
    2.5894448848121536e-05, // H+
    0.0062660162912796741, // HCO3-
    0.015741214773567921, // Ca+2
    0.0024602709470074127, // SO4-2
    1.8564927944498291,    // Cl-
    0.0099316034080546619, // Mg+2
    1.0725251492409775     // Na+1
  };

  timeStepTest< double, true >( carbonateSystem,
                                0.2,
                                10,
                                initialAggregateSpeciesConcentration,
                                expectedSpeciesConcentrations );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
