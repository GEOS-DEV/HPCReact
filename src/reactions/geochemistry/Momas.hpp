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

namespace geochemistry
{
// turn off uncrustify to allow for better readability of the parameters
// *****UNCRUSTIFY-OFF******

namespace momas
{

constexpr CArrayWrapper<double, 7, 12> stoichMatrix = 
  {
    // C1   C2   C3   C4   C5   CS1  CS2  |  X1  X2  X3   X4  S
    { -1,   0,   0,   0,   0,   0,   0,      0,  -1,  0,  0,  0 },  // C1 = -X2
    {  0,  -1,   0,   0,   0,   0,   0,      0,   1,  1,  0,  0 },  // C2 = X2 + X3
    {  0,   0,  -1,   0,   0,   0,   0,      0,  -1,  0,  1,  0 },  // C3 = X2 + X4
    {  0,   0,   0,  -1,   0,   0,   0,      0,  -4,  1,  3,  0 },  // C4 = -4X2 + X3 + 3X4
    {  0,   0,   0,   0,  -1,   0,   0,      0,   4,  3,  1,  0 },  // C5 = 4X2 + 3X3 + X4
    {  0,   0,   0,   0,   0,  -1,   0,      0,   3,  1,  0,  1 },  // CS1 = 3X2 + X3 + S
    {  0,   0,   0,   0,   0,   0,  -1,      0,  -3,  0,  1,  2 }   // CS2 = -3X2 + X4 + 2S
  };

constexpr CArrayWrapper<double, 7> equilibriumConstants = 
  { 
    1.0e12,   //   C1 + X2 = ??
    1.0,      //        C2 = X2 + X3
    1.0,      //        C3 = X2 + X4
    1.0e1,    //  C4 + 4X2 = X3 + 3X4
    1.0e-35,  //        C5 = 4X2 + 3X3 + X4
    1.0e-6,   //       CS1 = 3X2 + X3 + S
    1.0e1     // CS2 + 3X2 = + X4 + 2S
  };

constexpr CArrayWrapper<double, 7> forwardRates = 
  { 
    0.0,   // C1 = -X2
    0.0,   // C2 = X2 + X3
    0.0,   // C3 = X2 + X4
    0.0,   // C4 = -4X2 + X3 + 3X4
    0.0,   // C5 = 4X2 + 3X3 + X4
    0.0,   // CS1 = 3X2 + X3 + S
    0.0    // CS2 = -3X2 + X4 + 2S
  };

constexpr CArrayWrapper<double, 7> reverseRates = 
  { 
    0.0,   // C1 = -X2
    0.0,   // C2 = X2 + X3
    0.0,   // C3 = X2 + X4
    0.0,   // C4 = -4X2 + X3 + 3X4
    0.0,   // C5 = 4X2 + 3X3 + X4
    0.0,   // CS1 = 3X2 + X3 + S
    0.0    // CS2 = -3X2 + X4 + 2S
  };

constexpr CArrayWrapper<int, 7> mobileSpeciesFlag = 
  { 
    1,   // C1 = -X2
    1,   // C2 = X2 + X3
    1,   // C3 = X2 + X4
    1,   // C4 = -4X2 + X3 + 3X4
    1,   // C5 = 4X2 + 3X3 + X4
    0,   // CS1 = 3X2 + X3 + S
    0    // CS2 = -3X2 + X4 + 2S
  };
  }
  using momasSystemAllKineticType     = reactionsSystems::MixedReactionsParameters< double, int, int, 12, 7, 0 >;
  using momasSystemAllEquilibriumType = reactionsSystems::MixedReactionsParameters< double, int, int, 12, 7, 7 >;
  using momasSystemType               = reactionsSystems::MixedReactionsParameters< double, int, int, 12, 7, 7 >;

  constexpr momasSystemAllKineticType momasSystemAllKinetic( momas::stoichMatrix, momas::equilibriumConstants, momas::forwardRates, momas::reverseRates, momas::mobileSpeciesFlag );
  constexpr momasSystemAllEquilibriumType momasSystemAllEquilibrium( momas::stoichMatrix, momas::equilibriumConstants, momas::forwardRates, momas::reverseRates, momas::mobileSpeciesFlag );
  constexpr momasSystemType carbonateSystem( momas::stoichMatrix, momas::equilibriumConstants, momas::forwardRates, momas::reverseRates, momas::mobileSpeciesFlag );

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact
