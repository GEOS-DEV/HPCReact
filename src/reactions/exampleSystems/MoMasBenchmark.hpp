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
namespace MoMasBenchmark
{
// *****UNCRUSTIFY-OFF******

  using easyCaseType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 12, 7, 7 >;
  using mediumCaseType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 14, 10, 9 >; 

  constexpr easyCaseType easyCaseParams =
{
  // Stoichiometric matrix [7 rows × 12 columns]
  // Columns 0–6 are secondary species (must be -1 on the diagonal)
  // Columns 7–11 are primary species: {X1, X2, X3, X4, S}
  {
    // C1   C2   C3   C4   C5   CS1  CS2  |  X1  X2  X3   X4  S
    { -1,   0,   0,   0,   0,   0,   0,      0,  -1,  0,  0,  0 },  // C1 = -X2
    {  0,  -1,   0,   0,   0,   0,   0,      0,   1,  1,  0,  0 },  // C2 = X2 + X3
    {  0,   0,  -1,   0,   0,   0,   0,      0,  -1,  0,  1,  0 },  // C3 = -X2 + X4
    {  0,   0,   0,  -1,   0,   0,   0,      0,  -4,  1,  3,  0 },  // C4 = -4X2 + X3 + 3X4
    {  0,   0,   0,   0,  -1,   0,   0,      0,   4,  3,  1,  0 },  // C5 = 4X2 + 3X3 + X4
    {  0,   0,   0,   0,   0,  -1,   0,      0,   3,  1,  0,  1 },  // CS1 = 3X2 + X3 + S
    {  0,   0,   0,   0,   0,   0,  -1,      0,  -3,  0,  1,  2 }   // CS2 = -3X2 + X4 + 2S
  },

  // Equilibrium constants K
  {
    1.0e12,   //   C1 + X2 = inf
    1.0,      //        C2 = X2 + X3
    1.0,      //        C3 = -X2 + X4
    1.0e1,    //  C4 + 4X2 = X3 + 3X4
    1.0e-35,  //        C5 = 4X2 + 3X3 + X4
    1.0e-6,   //       CS1 = 3X2 + X3 + S
    1.0e1     // CS2 + 3X2 = + X4 + 2S
  },

  // Forward rate constants 
  { 
    0.0,   // C1 = -X2
    0.0,   // C2 = X2 + X3
    0.0,   // C3 = -X2 + X4
    0.0,   // C4 = -4X2 + X3 + 3X4
    0.0,   // C5 = 4X2 + 3X3 + X4
    0.0,   // CS1 = 3X2 + X3 + S
    0.0    // CS2 = -3X2 + X4 + 2S
  },

  // Reverse rate constants 
  { 
    0.0,   // C1 = -X2
    0.0,   // C2 = X2 + X3
    0.0,   // C3 = -X2 + X4
    0.0,   // C4 = -4X2 + X3 + 3X4
    0.0,   // C5 = 4X2 + 3X3 + X4
    0.0,   // CS1 = 3X2 + X3 + S
    0.0    // CS2 = -3X2 + X4 + 2S
  },

  // Flag of mobile secondary species
  { 1, // C1 = -X2
    1, // C2 = X2 + X3
    1, // C3 = -X2 + X4
    1, // C4 = -4X2 + X3 + 3X4
    1, // C5 = 4X2 + 3X3 + X4
    0, // CS1 = 3X2 + X3 + S
    0  // CS2 = -3X2 + X4 + 2S
  }
};

constexpr mediumCaseType mediumCaseParams =
{
  // Stoichiometric matrix [10 rows × 14 columns]
  // Columns 0–8 are secondary species (must be -1 on the diagonal)
  // Columns 9–13 are primary species: {X1, X2, X3, X4, S}
  {
    // C1   C2   C3   C4   C5   C6   C7   CS1  CS2  |  X1  X2  X3   X4  S
    { -1,   0,   0,   0,   0,   0,   0,   0,   0,      0,  -1,  0,  0,  0 },  // C1 = -X2
    {  0,  -1,   0,   0,   0,   0,   0,   0,   0,      0,   1,  1,  0,  0 },  // C2 = X2 + X3
    {  0,   0,  -1,   0,   0,   0,   0,   0,   0,      0,  -1,  0,  1,  0 },  // C3 = -X2 + X4
    {  0,   0,   0,  -1,   0,   0,   0,   0,   0,      0,  -4,  1,  3,  0 },  // C4 = -4X2 + X3 + 3X4
    {  0,   0,   0,   0,  -1,   0,   0,   0,   0,      0,   4,  3,  1,  0 },  // C5 = 4X2 + 3X3 + X4
    {  0,   0,   0,   0,   0,  -1,   0,   0,   0,      0,  10,  3,  0,  0 },  // C6 = 10X2 + 3X3
    {  0,   0,   0,   0,   0,   0,  -1,   0,   0,      0,  -8,  0,  2,  0 },  // C7 = -8X2 + 2X4
    {  0,   0,   0,   0,   0,   0,   0,  -1,   0,      0,   3,  1,  0,  1 },  // CS1 = 3X2 + X3 + S
    {  0,   0,   0,   0,   0,   0,   0,   0,  -1,      0,  -3,  0,  1,  2 },  // CS2 = -3X2 + X4 + 2S
    {  0,   0,   0,   0,   0,   0,   0,   0,   0,      0,  -3,  0,  1,  0 },  // Cc = -3X2 + X4 (kinetic)
  },

  // Equilibrium constants K
  {
    1.0e12,   //   C1 + X2 = inf
    1.0,      //        C2 = X2 + X3
    1.0,      //        C3 = -X2 + X4
    1.0e1,    //  C4 + 4X2 = X3 + 3X4
    1.0e-35,  //        C5 = 4X2 + 3X3 + X4
    1.0e-32,  //        C6 = 10X2 + 3X3
    1.0e4,    //  C7 + 8X2 = 2X4
    1.0e-6,   //       CS1 = 3X2 + X3 + S
    1.0e1,    // CS2 + 3X2 = X4 + 2S
    5         //  Cc + 3X2 = X4 (kinetic)
  },

  // Forward rate constants 
  { 
    0.0,   // C1 = -X2
    0.0,   // C2 = X2 + X3
    0.0,   // C3 = -X2 + X4
    0.0,   // C4 = -4X2 + X3 + 3X4
    0.0,   // C5 = 4X2 + 3X3 + X4
    0.0,   // C6 = 10X2 + 3X3
    0.0,   // C7 = -8X2 + 2X4
    0.0,   // CS1 = 3X2 + X3 + S
    0.0,   // CS2 = -3X2 + X4 + 2S
    10.0    // Cc = -3X2 + X4 (kinetic)
  },

  // Reverse rate constants 
  { 
    0.0,   // C1 = -X2
    0.0,   // C2 = X2 + X3
    0.0,   // C3 = -X2 + X4
    0.0,   // C4 = -4X2 + X3 + 3X4
    0.0,   // C5 = 4X2 + 3X3 + X4
    0.0,   // C6 = 10X2 + 3X3
    0.0,   // C7 = -8X2 + 2X4
    0.0,   // CS1 = 3X2 + X3 + S
    0.0,   // CS2 = -3X2 + X4 + 2S
    2.0    // Cc = -3X2 + X4 (kinetic)
  },

  // Flag of mobile secondary species
  { 1, // C1 = -X2
    1, // C2 = X2 + X3
    1, // C3 = -X2 + X4
    1, // C4 = -4X2 + X3 + 3X4
    1, // C5 = 4X2 + 3X3 + X4
    1, // C6 = 10X2 + 3X3
    1, // C7 = -8X2 + 2X4
    0, // CS1 = 3X2 + X3 + S
    0, // CS2 = -3X2 + X4 + 2S
    1  // Cc = -3X2 + X4 (kinetic)
  }
};

// *****UNCRUSTIFY-ON******
} // namespace MoMasBenchmark
} // namespace hpcReact
