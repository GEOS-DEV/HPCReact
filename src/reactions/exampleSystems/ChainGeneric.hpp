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

  using serialAllKineticType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 3, 3, 0 >;

  constexpr CArrayWrapper< signed char, 3, 3 > serialAllKineticStoichMatrix =
  {
    // Stoichiometric matrix [3 rows × 3 columns]
    // Columns 0–3 are primary species: {C1, C2, C3 }
    {
      // C1   C2   C3
      {  -1,   1,   0 },  // C1 = C2
      {   0,  -1,   1 },  // C2 = C3
      {   0,   0,  -1 },  // C3 =
    }
  };

  constexpr CArrayWrapper< double, 3 > serialAllKineticEquilibriumConstants =
  {
    0, // C1 = C2
    0, // C2 = C3
    0  // C3
  };

  constexpr CArrayWrapper< double, 3 > serialAllKineticForwardRates =
  {
    0.05, // C1 = C2
    0.03, // C2 = C3
    0.02  // C3
  };

  constexpr CArrayWrapper< double, 3 > serialAllKineticReverseRates =
  {
    0.0, // C1 = C2
    0.0, // C2 = C3
    0.0  // C3
  };

  constexpr CArrayWrapper< int, 3 > serialAllKineticMobileSpeciesFlag =
  {
    1, // C1 = C2
    1, // C2 = C3
    1  // C3
  };

  constexpr serialAllKineticType serialAllKineticParams(
    serialAllKineticStoichMatrix,
    serialAllKineticEquilibriumConstants,
    serialAllKineticForwardRates,
    serialAllKineticReverseRates,
    serialAllKineticMobileSpeciesFlag,
    0 );

  using serialAllKineticIonicStrengthType = SpeciatedIonicStrength< double, int, 3 >;

  using serialAllKineticActivityParamsType = Bdot< double, int, serialAllKineticIonicStrengthType >::Params;

  constexpr CArrayWrapper< double, 3 > serialAllKineticSpeciesCharge =
  { 0.0, 0.0, 0.0 };

  constexpr CArrayWrapper< double, 3 > serialAllKineticIonSize =
  { 3.5, 3.5, 3.5 };

  constexpr CArrayWrapper< double, 3 > serialAllKineticBdotParameters =
  { 0.0, 0.0, 0.0 };

  constexpr serialAllKineticActivityParamsType serialAllKineticActivityParams =
  {
    // species charge
    {{ serialAllKineticSpeciesCharge }},
    // ion size parameter
    serialAllKineticIonSize,
    // bdot parameter
    serialAllKineticBdotParameters
  };

  Identity< double, int, serialAllKineticIonicStrengthType >::Params serialAllKineticIdentityActivityParams = {};


// *****UNCRUSTIFY-ON******
} // namespace ChainGeneric
} // namespace hpcReact
