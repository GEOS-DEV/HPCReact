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

#include "reactions/unitTestUtilities/mixedReactionsTestUtilities.hpp"
#include "../GeochemicalSystems.hpp"
#include "constitutive/activity/Bdot.hpp"
#include "constitutive/activity/Identity.hpp"
#include "constitutive/ionicStrength/SpeciatedIonicStrength.hpp"


using namespace hpcReact;
using namespace hpcReact::unitTest_utilities;


/**
 * @brief Run the carbonate time step for a given activity model.
 * @details The system and initial state are identical across activity models, so only the model,
 *          its parameters, and the expected result vary between the tests below.
 */
template< typename ACTIVITY_MODEL >
void timeStepCarbonateSystemHelper( typename ACTIVITY_MODEL::Params const & activityParams,
                                    double const (&expectedSpeciesConcentrations)[hpcReact::geochemistry::carbonateSystemType::numPrimarySpecies()],
                                    double const calciteSurfaceArea,
                                    double const relativeTolerance = 1.0e-8 )
{
  using namespace hpcReact::geochemistry;

  static constexpr int numPrimarySpecies = carbonateSystemType::numPrimarySpecies();

  double const surfaceArea[carbonateSystemType::numKineticReactions()] =
  {
    calciteSurfaceArea, // CaCO3
  };

  double const initialAggregateSpeciesConcentration[numPrimarySpecies] =
  {
    3.76e-1, // H+
    3.76e-1, // HCO3-
    3.87e-2, // Ca+2
    3.21e-2, // SO4-2
    1.89, // Cl-
    1.65e-2, // Mg+2
    1.09 // Na+1
  };

  timeStepTest< double, true, ACTIVITY_MODEL >( carbonateSystem,
                                                activityParams,
                                                1.0,
                                                10,
                                                initialAggregateSpeciesConcentration,
                                                surfaceArea,
                                                expectedSpeciesConcentrations,
                                                relativeTolerance );
}


TEST( testMixedReactions, testTimeStep_carbonateSystem_Identity )
{
  using namespace hpcReact::geochemistry;

  // The Identity model leaves activities equal to concentrations, so these are the
  // ideal-solution concentrations.
  double const expectedSpeciesConcentrations[carbonateSystemType::numPrimarySpecies()] =
  {
    0.00043107371205575743, // H+
    0.00039393087146083564, // HCO3-
    0.015533203748051142, // Ca+2
    0.0025349216394956169, // SO4-2
    1.8597651352187075, // Cl-
    0.0099750254121169241, // Mg+2
    1.0720453048327119 // Na+1
  };

  timeStepCarbonateSystemHelper< carbonateNosolidIdentityActivityType >( carbonateNosolidIdentityActivityParams,
                                                                         expectedSpeciesConcentrations,
                                                                         1.0e-4 );
}


//******************************************************************************
// The B-dot activity model, verified against EQ6.
//
// EQ3NR speciates the brine and EQ6 then dissolves calcite into it under the same TST rate law,
// r = k*A*(1 - Q/K), so the expected values are its molalities after 10 s, read from
// eq36Database/calcite.6o. See the README there.
//
// The parameters are carbonateNosolidActivityParamsEQ36 rather than carbonateNosolidActivityParams:
// the comparison is only meaningful with the EQ3/6 B-dot parameters, since the phreeqc set leaves
// several of these complexes without an ion size and so at gamma = 1.
//
// The surface area is 0.01 m2 in HPCReact, which is the 100 cm2 of the EQ6 run.
//
// The tolerance is set by the reference, not by the solver, and matches the equilibrium test for
// the same reason: EQ6 prints five significant figures, and it fits the Debye-Huckel A and B over
// its temperature grid rather than reading the 25 C entry. Every primary species agrees to 1.4e-4.
TEST( testMixedReactions, testTimeStep_carbonateSystem_Bdot )
{
  using namespace hpcReact::geochemistry;

  double const expectedSpeciesConcentrations[carbonateSystemType::numPrimarySpecies()] =
  {
    1.3424e-04, // H+
    2.9929e-03, // HCO3-
    3.3722e-02, // Ca+2
    1.4972e-02, // SO4-2
    1.8834e+00, // Cl-
    1.4439e-02, // Mg+2
    1.0767e+00 // Na+1
  };

  timeStepCarbonateSystemHelper< carbonateNosolidActivityType >( carbonateNosolidActivityParamsEQ36,
                                                                 expectedSpeciesConcentrations,
                                                                 1.0e-2,
                                                                 5.0e-4 );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
