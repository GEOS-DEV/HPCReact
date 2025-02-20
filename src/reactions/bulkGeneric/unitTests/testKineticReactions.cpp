
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

// TEST( bulkGeneric, test_computeReactionRatesIntegral )
// {
//   using KineticReactionsType = KineticReactions< double, 
//                                                  double * const,
//                                                  double const * const, 
//                                                  int, 
//                                                  int >;

//   double const temperature = 298.15;
//   double primarySpeciesConcentration[] = { 1.0, 0.0, 0.5, 1.0, 0.0 };
//   double reactionRates[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

//   constexpr double dt = 1.0e-1;
//   constexpr int numPrimarySpecies = decltype(simpleTest)::numPrimarySpecies;
// //  constexpr int numKineticReactions = decltype(simpleTest)::numKineticReactions;

//   //double time = 0.0;
//   for( int i = 0; i < 1000; ++i )
//   {
//     KineticReactionsType::computeReactionRates( temperature, 
//                                                 simpleTest, 
//                                                 primarySpeciesConcentration, 
//                                                 reactionRates );

//     for( int j = 0; j < numPrimarySpecies; ++j )
//     {
//       primarySpeciesConcentration[j] += dt * reactionRates[j];
//     }
//     //time += dt;
//     // if( i % 10 == 0 )
//     // {
//     //   printf( " {% 12.8e} { % 12.8e, % 12.8e, % 12.8e, % 12.8e, % 12.8e }\n", 
//     //           time,
//     //           primarySpeciesConcentration[0],
//     //           primarySpeciesConcentration[1],
//     //           primarySpeciesConcentration[2],
//     //           primarySpeciesConcentration[3],
//     //           primarySpeciesConcentration[4] );
//     // }
//   }

//   EXPECT_NEAR( primarySpeciesConcentration[0], 3.92138294e-01, 1.0e-8 );
//   EXPECT_NEAR( primarySpeciesConcentration[1], 3.03930853e-01, 1.0e-8 );
//   EXPECT_NEAR( primarySpeciesConcentration[2], 5.05945481e-01, 1.0e-8 );
//   EXPECT_NEAR( primarySpeciesConcentration[3], 7.02014628e-01, 1.0e-8 );
//   EXPECT_NEAR( primarySpeciesConcentration[4], 5.95970745e-01, 1.0e-8 );

// }





TEST( bulkGeneric, test_computeReactionRates_ultramafic )
{
  using KineticReactionsType = KineticReactions< double, 
                                                 double * const,
                                                 double const * const, 
                                                 CArrayWrapper<double,5,5>,
                                                 int, 
                                                 int >;

  double const temperature = 298.15;
  constexpr int numPrimarySpecies = decltype(um1Params)::numSpecies;
  constexpr double minConc = 1.0e-6;
  //                                                        C02     H2CO3        H+     HCO3-    CO3^2-   Mg2SiO4      Mg2+      SiO2     MgCO3       H20
  double primarySpeciesConcentration[numPrimarySpecies] = { 10.0,     0.1,  minConc,  minConc, minConc,      20.0,  minConc,  minConc,  minConc,      1.0 };
  double reactionRates[numPrimarySpecies] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  CArrayWrapper<double,5,5> derivatives;
  constexpr double dt0 = 1.0e-14;

  double time = 0.0;
  for( int i = 0; i < 1; ++i )
  {
    KineticReactionsType::computeSpeciesRates( temperature, 
                                                um1Params, 
                                                primarySpeciesConcentration, 
                                                reactionRates,
                                                derivatives );


    double dt = dt0;
//     for( int j = 0; j < numPrimarySpecies; ++j )
//     {
//       double const delta = dt * reactionRates[j];
//       if( primarySpeciesConcentration[j] + delta < minConc )
//       {
//         double dtCut = ( - primarySpeciesConcentration[j] ) / reactionRates[j];
// //        printf( "reactionRates[%d], primarySpeciesConcentration[%d], dtCut = % 12.8e, % 12.8e, % 12.8e\n", j, j, reactionRates[j], primarySpeciesConcentration[j], dtCut );
//         if( dtCut < dt )
//         {
//           dt = dtCut;
//         } 
//       }
//     }
//    printf( "dt = % 12.8e\n", dt );

    for( int j = 0; j < numPrimarySpecies; ++j )
    {
      primarySpeciesConcentration[j] += dt * reactionRates[j];
    }
    time += dt;
      if( i == 0 )
      {         
        //       0         1         2         3         4         5         6         7         8         9
        //       01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
        //                 *         *         *         *         *         *         *         *         *         *         *
        printf( "      time        C02     H2CO3        H+     HCO3-    CO3^2-   Mg2SiO4      Mg2+      SiO2     MgCO3       H20\n");
      }
    if( i % 10000000 == 1 )
    {


        // printf( "           " );
        // for( int j = 0; j < numPrimarySpecies; ++j )
        // {
        //   printf( " % 6.2e", reactionRates[j] );
        // }
        // printf( "\n" );



        printf( " % 6.2e ", time );
        for( int j = 0; j < numPrimarySpecies; ++j )
        {
          printf( " % 6.2e", primarySpeciesConcentration[j] );
        }
        printf( "\n" );

    }
  }



}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
