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

#pragma once

#include "../reactionsSystems/Parameters.hpp"
#include "constitutive/ionicStrength/SpeciatedIonicStrength.hpp"
#include "constitutive/activity/Bdot.hpp"
#include "constitutive/activity/Identity.hpp"

namespace hpcReact
{
namespace ChainGeneric
{
// *****UNCRUSTIFY-OFF******

  using serialAllKineticType = reactionsSystems::MixedReactionsParameters< double, int, int, 3, 3, 0 >;

  constexpr serialAllKineticType serialAllKineticParams =
  {
    // Stoichiometric matrix [3 rows × 3 columns]
    // Columns 0–3 are primary species: {C1, C2, C3 }
    {
      // C1   C2   C3 
      {  -1,   1,   0 },  // C1 = C2
      {   0,  -1,   1 },  // C2 = C3
      {   0,   0,  -1 },  // C3 = 
    },

    // Equilibrium constants K
    {
      0, // C1 = C2
      0, // C2 = C3
      0  // C3
    },

    // Forward rate constants 
    { 
      0.05, // C1 = C2
      0.03, // C2 = C3
      0.02, // C3
    },

    // Reverse rate constants 
    { 
      0.0, // C1 = C2
      0.0, // C2 = C3 
      0.0  // C3
    },

    // Flag of mobile secondary species
    { 
      1, // C1 = C2
      1, // C2 = C3 
      1  // C3
    },

    0 // Use the forward and reverse to calculate the kinetic reaction rates
  };

  using serialAllKineticIonicStrengthType = SpeciatedIonicStrength< double, int, 3 >;

  using serialAllKineticActivityParamsType = Bdot< double, int, serialAllKineticIonicStrengthType >::Params;

  constexpr 
  serialAllKineticActivityParamsType
  serialAllKineticActivityParams =
  {
    // species charge
    {{ 0.0, 0.0, 0.0 }},
    // ion size parameter
    { 3.5, 3.5, 3.5 },
    // bdot parameter
    { 0.0, 0.0, 0.0 }
  };

  Identity< double, int, serialAllKineticIonicStrengthType >::Params serialAllKineticIdentityActivityParams = {};


// *****UNCRUSTIFY-ON******
} // namespace ChainGeneric
} // namespace hpcReact
