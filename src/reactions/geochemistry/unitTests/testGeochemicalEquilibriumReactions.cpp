
#include "reactions/unitTestUtilities/equilibriumReactionsTestUtilities.hpp"
#include "../GeochemicalSystems.hpp"

using namespace hpcReact;
using namespace hpcReact::geochemistry;
using namespace hpcReact::unitTest_utilities;


// //******************************************************************************
// TEST( testEquilibriumReactions, testEnforceEquilibrium )
// {
//   double const initialSpeciesConcentration[] = { 1.0, 1.0e-16, 0.5, 1.0, 1.0e-16 };
//   double const expectedSpeciesConcentrations[5] = { 3.92138294e-01, 3.03930853e-01, 5.05945481e-01, 7.02014628e-01, 5.95970745e-01 };


//   std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
//   testEnforceEquilibrium< double, 2 >( simpleTestRateParams.equilibriumReactionsParameters(),
//                                        initialSpeciesConcentration,
//                                        expectedSpeciesConcentrations );

// }


//******************************************************************************
TEST( testEquilibriumReactions, testcarbonateSystemAllEquilibrium )
{
  using namespace hpcReact::geochemistry;

  double const initialSpeciesConcentration[18] =
  {
    1.0e-16, // OH-
    1.0e-16, // CO2
    1.0e-16, // CO3-2
    1.0e-16, // H2CO3
    1.0e-16, // CaHCO3+
    1.0e-16, // CaSO4
    1.0e-16, // CaCl+
    1.0e-16, // CaCl2
    1.0e-16, // MgSO4
    1.0e-16, // NaSO4-
    1.0e-16, // CaCO3
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
    5.298329882666738e-03, // CaSO4
    5.844517547638333e-03, // CaCl+
    1.277319392670652e-02, // CaCl2
    6.618125707964991e-03, // MgSO4
    1.769217213462983e-02, // NaSO4-
    1.065032288527957e-09, // CaCO3
    4.396954721488358e-04, // H+
    3.723009698453808e-04, // HCO3-
    1.471656530812871e-02, // Ca+2
    2.491372274738741e-03, // SO4-2
    1.858609094598949e+00, // Cl-
    9.881874292035110e-03, // Mg+2
    1.072307827865370e+00 // Na+1
  };

  std::cout<<" RESIDUAL_FORM 0:"<<std::endl;
  testEnforceEquilibrium< double, 0 >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                       initialSpeciesConcentration,
                                       expectedSpeciesConcentrations );

  // std::cout<<" RESIDUAL_FORM 1:"<<std::endl;
  // testEnforceEquilibrium< double, 1 >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
  //                                      initialSpeciesConcentration,
  //                                      expectedSpeciesConcentrations );

  std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
  testEnforceEquilibrium< double, 2 >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                       initialSpeciesConcentration,
                                       expectedSpeciesConcentrations );

}


TEST( testEquilibriumReactions, testcarbonateSystemAllEquilibrium2 )
{

  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< double,
                                                                           int,
                                                                           int >;

  constexpr int numPrimarySpecies = hpcReact::geochemistry::carbonateSystemAllEquilibrium.numPrimarySpecies();

  double const initialPrimarySpeciesConcentration[numPrimarySpecies] =
  {
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89000, // Cl-
    1.65e-2, // Mg+2
    1.09000 // Na+1
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

  double logPrimarySpeciesConcentration[numPrimarySpecies];
  EquilibriumReactionsType::enforceEquilibrium_Aggregate( 0,
                                                          hpcReact::geochemistry::carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
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
