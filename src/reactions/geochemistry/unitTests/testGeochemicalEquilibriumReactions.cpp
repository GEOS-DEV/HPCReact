/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: (BSD-3-Clause)
 *
 * Copyright (c) 2025- Lawrence Livermore National Security LLC
 * All rights reserved
 *
 * See top level LICENSE files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


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
TEST( testEquilibriumReactions, testcarbonateSystemAllEquilibrium_Identity )
{
  using namespace hpcReact::geochemistry;

  double const initialSpeciesConcentration[17] =
  {
    1.0e-16, // OH-
    1.0e-16, // CO2
    1.0e-16, // CO3-2
    //1.0e-16, // H2CO3
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

  double const expectedSpeciesConcentrations[17] =
  { 1.7631991300262666e-11, // OH-
    0.37553965856049809, // CO2
    2.4208013881700254e-11, // CO3-2
    5.1052345859666575e-05, // CaHCO3+
    0.0050728241823463005, // CaSO4
    0.0058054702754571424, // CaCl+
    0.012170717821098593, // CaCl2
    0.0065270307598572714, // MgSO4
    0.017963901804350708, // NaSO4-
    0.00011324449476863112, // CaCO3
    0.00057358597611036207, // H+
    0.00029604457466600952, // HCO3-
    0.015486690880470161, // Ca+2
    0.0025362432534460182, // SO4-2
    1.8598530940823459, // Cl-
    0.0099729692401428292, // Mg+2
    1.0720360981956494 // Na+1
  };

  std::cout<<" RESIDUAL_FORM 0:"<<std::endl;
  testEnforceEquilibrium< double, 0, carbonateIdentityActivityType >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                                      hpcReact::geochemistry::carbonateIdentityActivityParams,
                                                                      initialSpeciesConcentration,
                                                                      expectedSpeciesConcentrations );

  // std::cout<<" RESIDUAL_FORM 1:"<<std::endl;
  // testEnforceEquilibrium< double, 1 >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
  //                                      initialSpeciesConcentration,
  //                                      expectedSpeciesConcentrations );

  std::cout<<" RESIDUAL_FORM 2:"<<std::endl;
  testEnforceEquilibrium< double, 2, carbonateIdentityActivityType >( carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                                      hpcReact::geochemistry::carbonateIdentityActivityParams,
                                                                      initialSpeciesConcentration,
                                                                      expectedSpeciesConcentrations );

}


TEST( testEquilibriumReactions, testcarbonateSystemAllEquilibrium2_Identity )
{


  static constexpr int numPrimarySpecies = hpcReact::geochemistry::carbonateSystemAllEquilibrium.numPrimarySpecies();
  static constexpr int numSpecies = hpcReact::geochemistry::carbonateSystemAllEquilibrium.numSpecies();

  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< double,
                                                                           int,
                                                                           int,
                                                                           Identity< double, int, SpeciatedIonicStrength< double, int, numSpecies > > >;

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
  EquilibriumReactionsType::enforceEquilibrium_PrimaryConcentrations( 0,
                                                                      hpcReact::geochemistry::carbonateSystemAllEquilibrium.equilibriumReactionsParameters(),
                                                                      hpcReact::geochemistry::carbonateIdentityActivityParams,
                                                                      initialPrimarySpeciesConcentration,
                                                                      logInitialPrimarySpeciesConcentration,
                                                                      logPrimarySpeciesConcentration );

  double const expectedPrimarySpeciesConcentrations[numPrimarySpecies] =
  {
    0.00057358597611063442, // H+
    0.00029604457466591314, // HCO3-
    0.015486690880470029, // Ca+2
    0.0025362432534459978, // SO4-2
    1.8598530940823459, // Cl-
    0.0099729692401427945, // Mg+2
    1.0720360981956494 // Na+1
  };

  for( int r=0; r<numPrimarySpecies; ++r )
  {
    EXPECT_NEAR( exp( logPrimarySpeciesConcentration[r] ), expectedPrimarySpeciesConcentrations[r], 1.0e-8 * expectedPrimarySpeciesConcentrations[r] );
  }


}

//******************************************************************************
// The B-dot activity model, verified against EQ3NR.
//
// The nine complexation reactions solved here are the model EQ3NR solves, so the expected values
// are its converged molalities, read from eq36Database/carbonate.3o. See the README there.
//
// The parameters are carbonateNosolidActivityParamsEQ36 rather than carbonateNosolidActivityParams:
// the comparison is only meaningful with the EQ3/6 B-dot parameters, since the phreeqc set leaves
// several of these complexes without an ion size and so at gamma = 1.
//
// The tolerance is set by the reference, not by the solver. EQ3NR prints five significant figures,
// and it fits the Debye-Huckel A and B over its temperature grid rather than reading the 25 C
// entry, so its effective values differ from the patched ones in the fifth decimal. Every primary
// species here agrees to 8.9e-5 or better.
TEST( testEquilibriumReactions, testcarbonateSystem_Bdot )
{
  using namespace hpcReact::geochemistry;

  static constexpr int numPrimarySpecies = carbonateSystemType::numPrimarySpecies();

  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< double, int, int,
                                                                           carbonateNosolidActivityType >;

  double const aggregatePrimarySpeciesConcentration[numPrimarySpecies] =
  {
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89, // Cl-
    1.65e-2, // Mg+2
    1.09 // Na+1
  };

  // The totals themselves, which is what a caller has before a run. B-dot cannot be started from
  // them directly; enforceEquilibrium_PrimaryConcentrations seeds itself with an ideal solve to get there.
  double logInitialGuess[numPrimarySpecies];
  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    logInitialGuess[i] = log( aggregatePrimarySpeciesConcentration[i] );
  }

  double logPrimarySpeciesConcentration[numPrimarySpecies];
  EquilibriumReactionsType::enforceEquilibrium_PrimaryConcentrations( 298.15,
                                                                      carbonateSystem.equilibriumReactionsParameters(),
                                                                      carbonateNosolidActivityParamsEQ36,
                                                                      aggregatePrimarySpeciesConcentration,
                                                                      logInitialGuess,
                                                                      logPrimarySpeciesConcentration );

  // EQ3NR converged molalities. This solve reports only the primary species; the secondary species
  // of the same run are compared in testCarbonateActivityVsEQ36.
  double const expectedPrimarySpeciesConcentrations[numPrimarySpecies] =
  {
    6.5867e-04, // H+
    6.1144e-04, // HCO3-
    3.2573e-02, // Ca+2
    1.4996e-02, // SO4-2
    1.8836e+00, // Cl-
    1.4435e-02, // Mg+2
    1.0766e+00 // Na+1
  };

  double const eq36Tolerance = 5.0e-4;

  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    EXPECT_NEAR( exp( logPrimarySpeciesConcentration[i] ),
                 expectedPrimarySpeciesConcentrations[i],
                 eq36Tolerance * expectedPrimarySpeciesConcentrations[i] );
  }
}

//******************************************************************************
// pH as the input, verified against EQ3NR.
//
// Same EQ3NR run as testcarbonateSystem_Bdot, read from eq36Database/carbonate.3o, but entered the
// way every reference code in the Xie suite enters it: pH in, total H out. The pH is the value on
// EQ3NR's "B-dot" scale, which is its internal scale (iopg(2) = -1, no rescaling), so it is
// -log10(a_H+) in exactly the sense this row constrains.
//
// The H+ total EQ3NR was given, 3.76e-1, is withheld here and must come back out. That closes the
// loop the other direction and is the check the mole balance rows cannot make.
//
// Every activity coefficient moves with the solve here, so this is the test that exercises the
// chain rule through the secondary species in the pX Jacobian row. It is also the only check on
// that row, and its tolerance is set by the reference rather than by the solver.
TEST( testEquilibriumReactions, testcarbonateSystem_Bdot_pHConstraint )
{
  using namespace hpcReact::geochemistry;
  using namespace hpcReact::reactionsSystems;

  static constexpr int numPrimarySpecies   = carbonateSystemType::numPrimarySpecies();
  static constexpr int numSecondarySpecies = carbonateSystemType::numSecondarySpecies();

  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< double, int, int,
                                                                           carbonateNosolidActivityType >;

  // eq36Database/carbonate.3o, "--- The pH, Eh, pe-, and Ah on various pH scales ---".
  constexpr double eq36pH = 3.2511;

  // The totals EQ3NR was given. Only H+ is replaced by the pH constraint; the rest stay as totals.
  double const aggregatePrimarySpeciesConcentration[numPrimarySpecies] =
  {
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89, // Cl-
    1.65e-2, // Mg+2
    1.09 // Na+1
  };

  double logInitialGuess[numPrimarySpecies];
  double constraintValue[numPrimarySpecies];
  PrimarySpeciesConstraintType constraintType[numPrimarySpecies];
  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    logInitialGuess[i] = log( aggregatePrimarySpeciesConcentration[i] );
    constraintValue[i] = aggregatePrimarySpeciesConcentration[i];
    constraintType[i] = PrimarySpeciesConstraintType::AggregateConcentration;
  }
  constraintValue[0] = eq36pH;
  constraintType[0] = PrimarySpeciesConstraintType::pX;

  double logPrimarySpeciesConcentration[numPrimarySpecies];
  double logSecondarySpeciesConcentration[numSecondarySpecies];
  double aggregatePrimarySpeciesConcentrationOut[numPrimarySpecies];
  EXPECT_TRUE( EquilibriumReactionsType::enforceEquilibrium_PrimaryConcentrations( 298.15,
                                                                                   carbonateSystem.equilibriumReactionsParameters(),
                                                                                   carbonateNosolidActivityParamsEQ36,
                                                                                   constraintType,
                                                                                   constraintValue,
                                                                                   logInitialGuess,
                                                                                   logPrimarySpeciesConcentration,
                                                                                   logSecondarySpeciesConcentration,
                                                                                   aggregatePrimarySpeciesConcentrationOut ) );

  // EQ3NR converged molalities, the same reference values testcarbonateSystem_Bdot compares against.
  double const expectedPrimarySpeciesConcentrations[numPrimarySpecies] =
  {
    6.5867e-04, // H+
    6.1144e-04, // HCO3-
    3.2573e-02, // Ca+2
    1.4996e-02, // SO4-2
    1.8836e+00, // Cl-
    1.4435e-02, // Mg+2
    1.0766e+00 // Na+1
  };

  // Set by the reference, not the solver: EQ3NR prints five figures and the pH is quoted to four.
  double const eq36Tolerance = 5.0e-4;

  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    EXPECT_NEAR( exp( logPrimarySpeciesConcentration[i] ),
                 expectedPrimarySpeciesConcentrations[i],
                 eq36Tolerance * expectedPrimarySpeciesConcentrations[i] );
  }

  EXPECT_NEAR( aggregatePrimarySpeciesConcentrationOut[0],
               aggregatePrimarySpeciesConcentration[0],
               eq36Tolerance * aggregatePrimarySpeciesConcentration[0] );
}

//******************************************************************************
// Charge balance on Cl-.
//
// When every other primary species is given its aggregate concentration T_p, electroneutrality
// fixes the Cl- aggregate concentration: T_Cl = -sum_{p != Cl} z_p T_p / z_Cl = 1.1362, where z_p
// is the charge. The test solves (a) with Cl- on charge balance and (b) with Cl- given 1.1362 as
// its aggregate concentration, and checks that the two results agree to roundoff.
TEST( testEquilibriumReactions, testcarbonateSystem_Bdot_chargeBalance )
{
  using namespace hpcReact::geochemistry;
  using namespace hpcReact::reactionsSystems;

  static constexpr int numPrimarySpecies   = carbonateSystemType::numPrimarySpecies();
  static constexpr int numSecondarySpecies = carbonateSystemType::numSecondarySpecies();

  using EquilibriumReactionsType = reactionsSystems::EquilibriumReactions< double, int, int,
                                                                           carbonateNosolidActivityType >;

  /// Index of Cl- among the primary species, the row that carries the balance.
  constexpr int balancedSpecies = 4;

  double const aggregatePrimarySpeciesConcentration[numPrimarySpecies] =
  {
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89, // Cl-
    1.65e-2, // Mg+2
    1.09 // Na+1
  };

  double logInitialGuess[numPrimarySpecies];
  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    logInitialGuess[i] = log( aggregatePrimarySpeciesConcentration[i] );
  }

  // The Cl- total that makes the water neutral, from the same identity the row uses.
  double chargeExcludingBalancedSpecies = 0.0;
  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    if( i != balancedSpecies )
    {
      chargeExcludingBalancedSpecies +=
        carbonateNosolidActivityParamsEQ36.m_speciesCharge[i + numSecondarySpecies] * aggregatePrimarySpeciesConcentration[i];
    }
  }
  double const neutralBalancedTotal =
    -chargeExcludingBalancedSpecies / carbonateNosolidActivityParamsEQ36.m_speciesCharge[balancedSpecies + numSecondarySpecies];

  EXPECT_NEAR( neutralBalancedTotal, 1.1362, 1.0e-12 );

  // (a) Cl- constrained by electroneutrality.
  PrimarySpeciesConstraintType constraintType[numPrimarySpecies];
  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    constraintType[i] = PrimarySpeciesConstraintType::AggregateConcentration;
  }
  constraintType[balancedSpecies] = PrimarySpeciesConstraintType::ChargeBalance;

  double logPrimaryBalanced[numPrimarySpecies];
  double logSecondaryBalanced[numSecondarySpecies];
  double aggregateBalanced[numPrimarySpecies];
  EXPECT_TRUE( EquilibriumReactionsType::enforceEquilibrium_PrimaryConcentrations( 298.15,
                                                                                   carbonateSystem.equilibriumReactionsParameters(),
                                                                                   carbonateNosolidActivityParamsEQ36,
                                                                                   constraintType,
                                                                                   aggregatePrimarySpeciesConcentration,
                                                                                   logInitialGuess,
                                                                                   logPrimaryBalanced,
                                                                                   logSecondaryBalanced,
                                                                                   aggregateBalanced ) );

  // (b) the same system with that total supplied directly.
  double constraintValue[numPrimarySpecies];
  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    constraintValue[i] = aggregatePrimarySpeciesConcentration[i];
    constraintType[i] = PrimarySpeciesConstraintType::AggregateConcentration;
  }
  constraintValue[balancedSpecies] = neutralBalancedTotal;

  double logPrimaryFromTotal[numPrimarySpecies];
  double logSecondaryFromTotal[numSecondarySpecies];
  double aggregateFromTotal[numPrimarySpecies];
  EXPECT_TRUE( EquilibriumReactionsType::enforceEquilibrium_PrimaryConcentrations( 298.15,
                                                                                   carbonateSystem.equilibriumReactionsParameters(),
                                                                                   carbonateNosolidActivityParamsEQ36,
                                                                                   constraintType,
                                                                                   constraintValue,
                                                                                   logInitialGuess,
                                                                                   logPrimaryFromTotal,
                                                                                   logSecondaryFromTotal,
                                                                                   aggregateFromTotal ) );

  for( int i = 0; i < numPrimarySpecies; ++i )
  {
    EXPECT_NEAR( exp( logPrimaryBalanced[i] ),
                 exp( logPrimaryFromTotal[i] ),
                 1.0e-12 * exp( logPrimaryFromTotal[i] ) );
  }

  EXPECT_NEAR( aggregateBalanced[balancedSpecies],
               neutralBalancedTotal,
               1.0e-10 * neutralBalancedTotal );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
