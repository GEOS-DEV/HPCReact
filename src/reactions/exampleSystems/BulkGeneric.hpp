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
namespace bulkGeneric
{
// turn off uncrustify to allow for better readability of the parameters
// *****UNCRUSTIFY-OFF******

// constexpr 
// EquilibriumKineticsModelConstants< double, int, int, 5 > 
// um1Constants = 
// {
//   { 1.70e-3, 3.2e-4, 4.7e-11,   1.0e7, 1.0e4 },
//   {   0.039,  1.0e4,   1.0e4, 1.0e-12, 1.0e8 },
//   {     0.0,    0.0,     0.0,     0.0,   0.0 },
// };

// constexpr 
// ReactionsParameters< double, int, int, EquilibriumKineticsModelConstants, 10, 5 > 
// um1Params = 
// { 
//   { 
//     { 
//       { -1,     1,  0,     0,      0,       0,    0,    0,     0,  -1 },
//       {  0,    -1,  1,     1,      0,       0,    0,    0,     0,   0 },
//       {  0,     0,  1,    -1,      1,       0,    0,    0,     0,   0 },
//       {  0,     0, -4,     0,      0,      -1,    2,    1,     0,   2 },
//       {  0,     0,  0,     0,     -1,       0,   -1,    0,     1,   0 } 
//     }
//   }, 
//   um1Constants };


using simpleKineticParamsType = reactionsSystems::KineticReactionsParameters< double, int, signed char, 5, 2 >;

constexpr CArrayWrapper< signed char, 2, 5 > simpleKineticStoichMatrix =
{
  // stoichiometric matrix
  {
    { -2, 1,  1,  0, 0 },
    {  0, 0, -1, -1, 2 }
  }
};

constexpr CArrayWrapper< double, 2 > simpleKineticForwardRates =
{ 1.0, 0.5 };

constexpr CArrayWrapper< double, 2 > simpleKineticReverseRates =
{ 1.0, 0.5 };

constexpr CArrayWrapper< double, 2 > simpleKineticEquilibriumConstants =
{ 1.0, 1.0 };

constexpr simpleKineticParamsType simpleKineticTestRateParams(
  simpleKineticStoichMatrix,
  simpleKineticForwardRates,
  simpleKineticReverseRates,
  simpleKineticEquilibriumConstants,
  0 );

using simpleIonicStrengthType = SpeciatedIonicStrength< double, int, 5 >;

using simpleActivityParamsType = Bdot< double, int, simpleIonicStrengthType >::Params;

constexpr CArrayWrapper< double, 5 > simpleSpeciesCharge =
{ 2.0, -1.0, 1.0, 0.0, -1.0 };

constexpr CArrayWrapper< double, 5 > simpleIonSize =
{ 4.0, 3.5, 3.5, 3.5, 3.5 };

constexpr CArrayWrapper< double, 5 > simpleBdotParameters =
{ 0.1, 0.1, 0.1, 0.0, 0.1 };

constexpr simpleActivityParamsType simpleActivityTestParams =
{
  // species charge
  {{ simpleSpeciesCharge }},
  // ion size parameter
  simpleIonSize,
  // bdot parameter
  simpleBdotParameters
};

Identity< double, int, simpleIonicStrengthType >::Params simpleIdentityActivityTestParams = {};


using simpleMixedParamsType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 5, 2, 2 >;

constexpr CArrayWrapper< signed char, 2, 5 > simpleMixedStoichMatrix =
{
  // stoichiometric matrix
  {
    { -2, 1,  1,  0, 0 },
    {  0, 0, -1, -1, 2 }
  }
};

constexpr CArrayWrapper< double, 2 > simpleMixedEquilibriumConstants =
{ 1.0, 1.0 };

constexpr CArrayWrapper< double, 2 > simpleMixedForwardRates =
{ 1.0, 0.5 };

constexpr CArrayWrapper< double, 2 > simpleMixedReverseRates =
{ 1.0, 0.5 };

constexpr CArrayWrapper< int, 2 > simpleMixedMobileSpeciesFlag =
{ 1, 1 };

constexpr simpleMixedParamsType simpleTestRateParams(
  simpleMixedStoichMatrix,
  simpleMixedEquilibriumConstants,
  simpleMixedForwardRates,
  simpleMixedReverseRates,
  simpleMixedMobileSpeciesFlag,
  0 );

// *****UNCRUSTIFY-ON******
} // namespace bulkGeneric
} // namespace hpcReact
