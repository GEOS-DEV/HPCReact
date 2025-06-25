#include "reactions/unitTestUtilities/mixedReactionsTestUtilities.hpp"
#include "../GeochemicalSystems.hpp"


using namespace hpcReact;
using namespace hpcReact::unitTest_utilities;


TEST( testMixedReactions, testTimeStep_carbonateSystem )
{
  using namespace hpcReact::geochemistry;

  constexpr int numPrimarySpecies = carbonateSystemType::numPrimarySpecies();

  double const surfaceArea[carbonateSystemType::numKineticReactions()] =
  {
    1.0, // CaCO3
  };

  double const initialAggregateSpeciesConcentration[numPrimarySpecies] =
  {
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
    9.671777755634228e-06, // H+
    0.016494147899655441,  // HCO3-
    0.017327206801111415,  // Ca+2
    0.0024137247729776557, // SO4-2
    1.8532341292597552,    // Cl-
    0.010006970034001514,  // Mg+2
    1.0728505565167725     // Na+1
  };

  timeStepTest< double, true >( carbonateSystem,
                                1.0e-8,
                                10,
                                initialAggregateSpeciesConcentration,
                                surfaceArea,
                                expectedSpeciesConcentrations );

}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
