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


using simpleKineticTestType = reactionsSystems::MixedReactionsParameters< double, int, int, 5, 2, 0 >;

constexpr 
simpleKineticTestType
simpleKineticTestRateParams = 
{ 
  // stoichiometric matrix
  {
    { -2, 1,  1,  0, 0 },
    {  0, 0, -1, -1, 2 }
  },
  // equilibrium constants
  { 1.0, 1.0 },
  // forward rate constants
  { 1.0, 0.5 },
  // reverse rate constants
  { 1.0, 0.5 }
};

using simpleTestType = reactionsSystems::MixedReactionsParameters< double, int, int, 5, 2, 2 >;

constexpr 
simpleTestType
simpleTestRateParams = 
{ 
  // stoichiometric matrix
  {
    { -2, 1,  1,  0, 0 },
    {  0, 0, -1, -1, 2 }
  },
  // equilibrium constants
  { 1.0, 1.0 },
  // forward rate constants
  { 1.0, 0.5 },
  // reverse rate constants
  { 1.0, 0.5 }
};

// *****UNCRUSTIFY-ON******
} // namespace bulkGeneric
} // namespace hpcReact
