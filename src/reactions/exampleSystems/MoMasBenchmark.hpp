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

  using easyCaseType = reactionsSystems::MixedReactionsParameters< double, int, int, 12, 7, 7 >;
  using mediumCaseType = reactionsSystems::MixedReactionsParameters< double, int, int, 14, 10, 9 >; 

  constexpr CArrayWrapper< double, 7, 12 > easyCaseStoichMatrix =
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
    }
  };

  constexpr CArrayWrapper< double, 7 > easyCaseEquilibriumConstants =
  {
    1.0e12,   //   C1 + X2 = inf
    1.0,      //        C2 = X2 + X3
    1.0,      //        C3 = -X2 + X4
    1.0e1,    //  C4 + 4X2 = X3 + 3X4
    1.0e-35,  //        C5 = 4X2 + 3X3 + X4
    1.0e-6,   //       CS1 = 3X2 + X3 + S
    1.0e1     // CS2 + 3X2 = + X4 + 2S
  };

  constexpr CArrayWrapper< double, 7 > easyCaseForwardRates =
  {
    0.0,   // C1 = -X2
    0.0,   // C2 = X2 + X3
    0.0,   // C3 = -X2 + X4
    0.0,   // C4 = -4X2 + X3 + 3X4
    0.0,   // C5 = 4X2 + 3X3 + X4
    0.0,   // CS1 = 3X2 + X3 + S
    0.0    // CS2 = -3X2 + X4 + 2S
  };

  constexpr CArrayWrapper< double, 7 > easyCaseReverseRates =
  {
    0.0,   // C1 = -X2
    0.0,   // C2 = X2 + X3
    0.0,   // C3 = -X2 + X4
    0.0,   // C4 = -4X2 + X3 + 3X4
    0.0,   // C5 = 4X2 + 3X3 + X4
    0.0,   // CS1 = 3X2 + X3 + S
    0.0    // CS2 = -3X2 + X4 + 2S
  };

  constexpr CArrayWrapper< int, 7 > easyCaseMobileSpeciesFlag =
  {
    1, // C1 = -X2
    1, // C2 = X2 + X3
    1, // C3 = -X2 + X4
    1, // C4 = -4X2 + X3 + 3X4
    1, // C5 = 4X2 + 3X3 + X4
    0, // CS1 = 3X2 + X3 + S
    0  // CS2 = -3X2 + X4 + 2S
  };

  constexpr easyCaseType easyCaseParams(
    easyCaseStoichMatrix,
    easyCaseEquilibriumConstants,
    easyCaseForwardRates,
    easyCaseReverseRates,
    easyCaseMobileSpeciesFlag );

  constexpr CArrayWrapper< double, 10, 14 > mediumCaseStoichMatrix =
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
    }
  };

  constexpr CArrayWrapper< double, 10 > mediumCaseEquilibriumConstants =
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
  };

  constexpr CArrayWrapper< double, 10 > mediumCaseForwardRates =
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
    10.0   // Cc = -3X2 + X4 (kinetic)
  };

  constexpr CArrayWrapper< double, 10 > mediumCaseReverseRates =
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
  };

  constexpr CArrayWrapper< int, 10 > mediumCaseMobileSpeciesFlag =
  {
    1, // C1 = -X2
    1, // C2 = X2 + X3
    1, // C3 = -X2 + X4
    1, // C4 = -4X2 + X3 + 3X4
    1, // C5 = 4X2 + 3X3 + X4
    1, // C6 = 10X2 + 3X3
    1, // C7 = -8X2 + 2X4
    0, // CS1 = 3X2 + X3 + S
    0, // CS2 = -3X2 + X4 + 2S
    1  // Cc = -3X2 + X4 (kinetic)
  };

  constexpr mediumCaseType mediumCaseParams(
    mediumCaseStoichMatrix,
    mediumCaseEquilibriumConstants,
    mediumCaseForwardRates,
    mediumCaseReverseRates,
    mediumCaseMobileSpeciesFlag );

// *****UNCRUSTIFY-ON******
} // namespace MoMasBenchmark
} // namespace hpcReact
