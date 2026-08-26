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
                                    double const (&expectedSpeciesConcentrations)[hpcReact::geochemistry::carbonateSystemType::numPrimarySpecies()] )
{
  using namespace hpcReact::geochemistry;

  static constexpr int numPrimarySpecies = carbonateSystemType::numPrimarySpecies();

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
                                                expectedSpeciesConcentrations );
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
                                                                         expectedSpeciesConcentrations );
}


//******************************************************************************
// Placeholder for the B-dot carbonate model.
//
// This one will be checked against EQ3/6 for verification.
//
// The values below are deliberately zero: if this is enabled before the EQ3/6 run is done, it
// fails immediately rather than appearing to pass.
TEST( testMixedReactions, DISABLED_testTimeStep_carbonateSystem_Bdot )
{
  using namespace hpcReact::geochemistry;

  double const expectedSpeciesConcentrations[carbonateSystemType::numPrimarySpecies()] =
  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // TODO: from EQ3/6

  timeStepCarbonateSystemHelper< carbonateNosolidActivityType >( carbonateNosolidActivityParams,
                                                                 expectedSpeciesConcentrations );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
