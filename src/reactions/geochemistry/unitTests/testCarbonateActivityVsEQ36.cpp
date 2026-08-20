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

/**
 * @file testCarbonateActivityVsEQ36.cpp
 * @brief Validation of the activity model against EQ3NR (EQ3/6 v8.0a).
 *
 * Reference data is generated with the 'cmp' database (data0.com.V8.R6), whose adh/bdh entries are
 * patched to the A and B this code derives from physical constants so that the only remaining
 * difference is the model itself. Species outside the carbonate system's 17-species set are
 * suppressed in the EQ3NR input. The database, input and output files are in eq36Database/;
 * see the README there to regenerate these values.
 *
 * These tests use carbonateActivityParamsEQ36.
 */

#include "../Carbonate.hpp"
#include "common/constants.hpp"

#include <gtest/gtest.h>

using namespace hpcReact;
using namespace hpcReact::geochemistry;

namespace
{

constexpr int numSpecies = 17;

/// EQ3NR converged molalities, in this code's species order.
constexpr double eq36Molality[numSpecies] =
{
  2.6702e-11, // OH-
  3.7534e-01, // CO2(aq)
  2.3166e-10, // CO3-2
  4.7225e-05, // CaHCO3+
  1.6821e-03, // CaSO4(aq)
  2.3750e-03, // CaCl+
  2.0222e-03, // CaCl2(aq)
  2.0650e-03, // MgSO4(aq)
  1.3357e-02, // NaSO4-
  5.0225e-10, // CaCO3(aq)
  6.5867e-04, // H+
  6.1144e-04, // HCO3-
  3.2573e-02, // Ca+2
  1.4996e-02, // SO4-2
  1.8836e+00, // Cl-
  1.4435e-02, // Mg+2
  1.0766e+00 // Na+
};

/// EQ3NR log10(gamma), printed to four decimals.
constexpr double eq36Log10Gamma[numSpecies] =
{
  -0.1965, // OH-
  0.1555, // CO2(aq)   (Drummond salting-out, not reproduced here)
  -0.8321, // CO3-2
  -0.1760, // CaHCO3+
  0.0000, // CaSO4(aq)
  -0.1760, // CaCl+
  0.0000, // CaCl2(aq)
  0.0000, // MgSO4(aq)
  -0.1760, // NaSO4-
  0.0000, // CaCO3(aq)
  -0.0698, // H+
  -0.1760, // HCO3-
  -0.6717, // Ca+2
  -0.9023, // SO4-2
  -0.2208, // Cl-
  -0.5298, // Mg+2
  -0.1760 // Na+
};

constexpr double eq36IonicStrength = 1.6126;

/// EQ3/6 gives CO2(aq) a Drummond salting-out term; this model gives every neutral gamma = 1.
constexpr int indexCO2 = 1;

} // namespace


TEST( testCarbonateActivityVsEQ36, ionicStrength )
{
  double dIonicStrength_dConcentration[numSpecies];
  double const I = carbonateIonicStrengthType::calculate( carbonateActivityParamsEQ36,
                                                          eq36Molality,
                                                          dIonicStrength_dConcentration );
  EXPECT_NEAR( I, eq36IonicStrength, 1.0e-4 );
}


TEST( testCarbonateActivityVsEQ36, activityCoefficients )
{
  double logActivityCoefficients[numSpecies] = { 0.0 };
  double dLogActivityCoefficients_dConcentrations[numSpecies][numSpecies] = {{ 0.0 }};

  carbonateActivityType::calculateLogActivityCoefficients( carbonateActivityParamsEQ36,
                                                           eq36Molality,
                                                           logActivityCoefficients,
                                                           dLogActivityCoefficients_dConcentrations );

  for( int i = 0; i < numSpecies; ++i )
  {
    if( i == indexCO2 )
      continue;
    // EQ3NR truncates log10(gamma) to four decimals, biasing the reference low by up to 1e-4.
    EXPECT_NEAR( logActivityCoefficients[i] * constants::invln10, eq36Log10Gamma[i], 2.0e-4 )
      << "species index " << i;
  }
}


TEST( testCarbonateActivityVsEQ36, neutralSpeciesAreIdeal )
{
  double logActivityCoefficients[numSpecies] = { 0.0 };
  double dLogActivityCoefficients_dConcentrations[numSpecies][numSpecies] = {{ 0.0 }};

  carbonateActivityType::calculateLogActivityCoefficients( carbonateActivityParamsEQ36,
                                                           eq36Molality,
                                                           logActivityCoefficients,
                                                           dLogActivityCoefficients_dConcentrations );

  EXPECT_NEAR( logActivityCoefficients[indexCO2], 0.0, 1.0e-12 );
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
