
#include "../KineticReactions_impl.hpp"
#include "../ParametersPredefined.hpp"
#include "common/DirectSystemSolve.hpp"

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
//   double speciesConcentration[6] = { 0.01, 0.01, 0.01, 0.01, 0.01, 1.0 };
//   double speciesRates[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  
//   KineticReactionsType::computeReactionRates( temperature, 
//                                               bicarbonateBuffer, 
//                                               speciesConcentration, 
//                                               speciesRates );

//   printf( "speciesRates = {%8.4e, %8.4e, %8.4e, %8.4e, %8.4e, %8.4e }\n", 
//           speciesRates[0],
//           speciesRates[1],
//           speciesRates[2],
//           speciesRates[3],
//           speciesRates[4],
//           speciesRates[5] );
//   // EXPECT_NEAR( speciesRates[0], 0.0, 1.0e-8 );
//   // EXPECT_NEAR( speciesRates[1], 0.0, 1.0e-8 );
//   // EXPECT_NEAR( speciesRates[2], 0.0, 1.0e-8 );
//   // EXPECT_NEAR( speciesRates[3], 0.0, 1.0e-8 );
//   // EXPECT_NEAR( speciesRates[4], 0.0, 1.0e-8 );
  
// }

TEST( bulkGeneric, test_computeReactionRatesIntegral )
{
  using KineticReactionsType = KineticReactions< double, 
                                                 int, 
                                                 int >;

  double const temperature = 298.15;
  double speciesConcentration[] = { 1.0, 0.0, 0.5, 1.0, 0.0 };
  double speciesRates[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

  constexpr double dt = 2.0;
  constexpr int numSpecies = decltype(simpleTestRateParams)::numSpecies;
//  constexpr int numReactions = decltype(simpleTestRateParams)::numReactions;


  double time = 0.0;
  for( int t = 0; t < 10; ++t )
  {
    CArrayWrapper<double,5,5> speciesRateDerivatives;
    double speciesConcentration_n[] = { speciesConcentration[0], 
                                               speciesConcentration[1], 
                                               speciesConcentration[2], 
                                               speciesConcentration[3], 
                                               speciesConcentration[4] };



    double residualNorm = 0.0;
    for( int k=0; k<10; ++k )
    {
      // printf( "iteration %2d: \n", k );

      KineticReactionsType::computeSpeciesRates( temperature, 
                                                 simpleTestRateParams, 
                                                 speciesConcentration, 
                                                 speciesRates,
                                                 speciesRateDerivatives );

      double residual[numSpecies] = { 0.0, 0.0, 0.0, 0.0, 0.0 };
      double deltaPrimarySpeciesConcentration[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

      for( int i = 0; i < numSpecies; ++i )
      {
        residual[i] = -(speciesConcentration[i] - speciesConcentration_n[i] - dt * speciesRates[i]);
        for( int j = 0; j < numSpecies; ++j )
        {
          speciesRateDerivatives( i, j ) = - dt * speciesRateDerivatives( i, j );
        }
        speciesRateDerivatives( i, i ) += 1.0;
      }

      residualNorm = 0.0;
      for( int j = 0; j < numSpecies; ++j )
      {
        residualNorm += residual[j] * residual[j];
      }
      residualNorm = sqrt( residualNorm );
      // printf( "     residual = { % 8.4e, % 8.4e, % 8.4e, % 8.4e, % 8.4e } -> ||% 8.4e||\n", 
      //         residual[0],
      //         residual[1],
      //         residual[2],
      //         residual[3],
      //         residual[4],
      //         residualNorm );
      if( residualNorm < 1.0e-8 )
      {
        break;
      }
      // printf(    "dR/dC :\n" );
      // for( int i = 0; i < numSpecies; ++i )
      // {
      //   printf( "                | % 6.3g, % 6.3g, % 6.3g, % 6.3g, % 6.3g |\n", 
      //           speciesRateDerivatives( i, 0 ),
      //           speciesRateDerivatives( i, 1 ),
      //           speciesRateDerivatives( i, 2 ),
      //           speciesRateDerivatives( i, 3 ),
      //           speciesRateDerivatives( i, 4 ) );
      // }

      solveNxN_pivoted<double,numSpecies>( speciesRateDerivatives.data, residual, deltaPrimarySpeciesConcentration );

      // printf( "       deltaC = { % 8.4e, % 8.4e, % 8.4e, % 8.4e, % 8.4e }\n", 
      //   deltaPrimarySpeciesConcentration[0],
      //   deltaPrimarySpeciesConcentration[1],
      //   deltaPrimarySpeciesConcentration[2],
      //   deltaPrimarySpeciesConcentration[3],
      //   deltaPrimarySpeciesConcentration[4] );
      for( int i = 0; i < numSpecies; ++i )
      {
        speciesConcentration[i] += deltaPrimarySpeciesConcentration[i];
      }
      // printf( "        C_k+1 = { % 8.4e, % 8.4e, % 8.4e, % 8.4e, % 8.4e }\n\n", 
      //   speciesConcentration[0],
      //   speciesConcentration[1],
      //   speciesConcentration[2],
      //   speciesConcentration[3],
      //   speciesConcentration[4] );


    }

    time += dt;
//    if( i % 10 == 0 )
    {
      printf( " time : residual = %4.1f : % 8.4e { % 8.4e, % 8.4e, % 8.4e, % 8.4e, % 8.4e }\n", 
              time, residualNorm,
              speciesConcentration[0],
              speciesConcentration[1],
              speciesConcentration[2],
              speciesConcentration[3],
              speciesConcentration[4] );
    }
  }

  EXPECT_NEAR( speciesConcentration[0], 3.92138294e-01, 1.0e-4 );
  EXPECT_NEAR( speciesConcentration[1], 3.03930853e-01, 1.0e-4 );
  EXPECT_NEAR( speciesConcentration[2], 5.05945481e-01, 1.0e-4 );
  EXPECT_NEAR( speciesConcentration[3], 7.02014628e-01, 1.0e-4 );
  EXPECT_NEAR( speciesConcentration[4], 5.95970745e-01, 1.0e-4 );

}

TEST( bulkGeneric, test_timeStep )
{
  using KineticReactionsType = KineticReactions< double, 
                                                 int, 
                                                 int >;

  double const temperature = 298.15;
  double speciesConcentration[] = { 1.0, 0.0, 0.5, 1.0, 0.0 };
  double speciesRates[] = { 0.0, 0.0, 0.0, 0.0, 0.0 };

  constexpr double dt = 2.0;
  constexpr int numSpecies = decltype(simpleTestRateParams)::numSpecies;


  double time = 0.0;
  for( int t = 0; t < 10; ++t )
  {
    CArrayWrapper<double,numSpecies,numSpecies> speciesRateDerivatives;
    double speciesConcentration_n[] = { speciesConcentration[0], 
                                               speciesConcentration[1], 
                                               speciesConcentration[2], 
                                               speciesConcentration[3], 
                                               speciesConcentration[4] };



    KineticReactionsType::timeStep( dt, 
                                    temperature, 
                                    simpleTestRateParams, 
                                    speciesConcentration_n, 
                                    speciesConcentration, 
                                    speciesRates,
                                    speciesRateDerivatives );

    time += dt;
    
    {
      printf( " time = %4.1f { % 8.4e, % 8.4e, % 8.4e, % 8.4e, % 8.4e }\n", 
              time,
              speciesConcentration[0],
              speciesConcentration[1],
              speciesConcentration[2],
              speciesConcentration[3],
              speciesConcentration[4] );
    }
  }

  EXPECT_NEAR( speciesConcentration[0], 3.92138294e-01, 1.0e-4 );
  EXPECT_NEAR( speciesConcentration[1], 3.03930853e-01, 1.0e-4 );
  EXPECT_NEAR( speciesConcentration[2], 5.05945481e-01, 1.0e-4 );
  EXPECT_NEAR( speciesConcentration[3], 7.02014628e-01, 1.0e-4 );
  EXPECT_NEAR( speciesConcentration[4], 5.95970745e-01, 1.0e-4 );

}




// TEST( bulkGeneric, test_computeReactionRates_ultramafic )
// {
//   using KineticReactionsType = KineticReactions< double, 
//                                                  int, 
//                                                  int >;

//   double const temperature = 298.15;
//   constexpr int numPrimarySpecies = decltype(um1Params)::numSpecies;
//   constexpr double minConc = 1.0e-6;
//   //                                                        C02     H2CO3        H+     HCO3-    CO3^2-   Mg2SiO4      Mg2+      SiO2     MgCO3       H20
//   double speciesConcentration[numPrimarySpecies] = { 10.0,     0.1,  minConc,  minConc, minConc,      20.0,  minConc,  minConc,  minConc,      1.0 };
//   double speciesRates[numPrimarySpecies] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
//   CArrayWrapper<double,5,5> derivatives;
//   constexpr double dt0 = 1.0e-14;

//   double time = 0.0;
//   for( int i = 0; i < 1; ++i )
//   {
//     KineticReactionsType::computeSpeciesRates( temperature, 
//                                                 um1Params, 
//                                                 speciesConcentration, 
//                                                 speciesRates,
//                                                 derivatives );


//     double dt = dt0;
// //     for( int j = 0; j < numPrimarySpecies; ++j )
// //     {
// //       double const delta = dt * speciesRates[j];
// //       if( speciesConcentration[j] + delta < minConc )
// //       {
// //         double dtCut = ( - speciesConcentration[j] ) / speciesRates[j];
// // //        printf( "speciesRates[%d], speciesConcentration[%d], dtCut = % 8.4e, % 12.8e, % 12.8e\n", j, j, speciesRates[j], speciesConcentration[j], dtCut );
// //         if( dtCut < dt )
// //         {
// //           dt = dtCut;
// //         } 
// //       }
// //     }
// //    printf( "dt = % 12.8e\n", dt );

//     for( int j = 0; j < numPrimarySpecies; ++j )
//     {
//       speciesConcentration[j] += dt * speciesRates[j];
//     }
//     time += dt;
//       if( i == 0 )
//       {         
//         //       0         1         2         3         4         5         6         7         8         9
//         //       01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
//         //                 *         *         *         *         *         *         *         *         *         *         *
//         printf( "      time        C02     H2CO3        H+     HCO3-    CO3^2-   Mg2SiO4      Mg2+      SiO2     MgCO3       H20\n");
//       }
//     if( i % 10000000 == 1 )
//     {


//         // printf( "           " );
//         // for( int j = 0; j < numPrimarySpecies; ++j )
//         // {
//         //   printf( " % 6.2e", speciesRates[j] );
//         // }
//         // printf( "\n" );



//         printf( " % 6.2e ", time );
//         for( int j = 0; j < numPrimarySpecies; ++j )
//         {
//           printf( " % 6.2e", speciesConcentration[j] );
//         }
//         printf( "\n" );

//     }
//   }



// }

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
