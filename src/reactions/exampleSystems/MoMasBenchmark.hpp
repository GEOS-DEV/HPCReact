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
namespace MomMasBenchmark
{
// *****UNCRUSTIFY-OFF******

  using simpleSystemType = reactionsSystems::MixedReactionsParameters< double, int, int, 12, 7, 7 >;

  constexpr simpleSystemType simpleSystemParams =
{
  // Stoichiometric matrix [7 rows × 12 columns]
  // Columns 0–6 are secondary species (must be -1 on the diagonal)
  // Columns 7–11 are primary species: {X1, X2, X3, X4, S}
  {
    // C1   C2   C3   C4   C5   CS1  CS2  |  X1  X2  X3   X4  S
    { -1,   0,   0,   0,   0,   0,   0,      0,  -1,  0,  0,  0 },  // C1 = X2
    {  0,  -1,   0,   0,   0,   0,   0,      0,   1,  1,  0,  0 },  // C2 = X2 * X3
    {  0,   0,  -1,   0,   0,   0,   0,      0,  -1,  0,  1,  0 },  // C3 = X2 * X4
    {  0,   0,   0,  -1,   0,   0,   0,      0,  -4,  1,  3,  0 },  // C4 = -4X2 + X3 + 3X4
    {  0,   0,   0,   0,  -1,   0,   0,      0,   4,  3,  1,  0 },  // C5 = 4X2 + 3X3 + X4
    {  0,   0,   0,   0,   0,  -1,   0,      0,   3,  1,  0,  1 },  // CS1 = 3X2 + X3 + S
    {  0,   0,   0,   0,   0,   0,  -1,      0,  -3,  0,  1,  2 }   // CS2 = -3X2 + X4 + 2S
  },

  // Equilibrium constants K
  {
    1.0e-12,   // C1
    1.0,       // C2
    1.0,       // C3
    0.1,       // C4
    1.0e35,    // C5
    1.0e6,     // CS1
    1.0e-1     // CS2
  },

  // Forward rate constants 
  { 0.0,
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0 
  },

  // Reverse rate constants 
  { 0.0,
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0 
  },

  // Flag of mobile secondary species
  { 1,
    1,
    1,
    1,
    1,
    0,
    0
  }
};

// *****UNCRUSTIFY-ON******
} // namespace MomMasBenchmark
} // namespace hpcReact
