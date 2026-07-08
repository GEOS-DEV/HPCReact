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

namespace mascagnite
{

constexpr CArrayWrapper<signed char, 3, 5> stoichMatrix = 
  { //   NH3(aq) HSO4-  H+    NH4+  SO4--   
    {     1,     0,     1,   -1,    0 }, //  NH4+ = NH3(aq) + H+ 
    {     0,    -1,     1,    0 ,   1 }, //  HSO4- = H+ + SO4--
    {     0,     0,     0,    2,    1 } //  (NH4)2SO4(s) = 2NH4+ + SO4--
  };

  constexpr CArrayWrapper<signed char, 1, 2> stoichMatrixNoAqueous = 
  { // NH4+  SO4--   
    {  2,   1 } //  (NH4)2SO4(s) = 2NH4+ + SO4--
  };

constexpr CArrayWrapper<double, 3> equilibriumConstants = 
  { 
    5.59757601e-10,  //  NH4+ = NH3(aq) + H+ 
    0.01028016,      //  HSO4- = H+ + SO4--
    2.0328        //  (NH4)2SO4(s) = 2NH4+ + SO4--
  };

constexpr CArrayWrapper<double, 3> forwardRates = 
  { 
    1.0e10,  //  NH4+ = NH3(aq) + H+ 
    1.0e10,  //  HSO4- = H+ + SO4--
    0.3      //  (NH4)2SO4(s) = 2NH4+ + SO4--
  };

constexpr CArrayWrapper<double, 3> reverseRates = 
  { 
    1.0e10,      //  NH4+ = NH3(aq) + H+ 
    1.0e10,      //  HSO4- = H+ + SO4--
    0.1475   //  (NH4)2SO4(s) = 2NH4+ + SO4--
  };

constexpr CArrayWrapper<int, 3> mobileSpeciesFlag = 
  { 
    
    1,  //  NH4+ = NH3(aq) + H+ 
    1,  //  HSO4- = H+ + SO4--
    1   //  (NH4)2SO4(s) = 2NH4+ + SO4--
  };

}

using mascagniteSystemType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 5, 3, 2 >;

constexpr mascagniteSystemType mascagniteSystem( mascagnite::stoichMatrix, mascagnite::equilibriumConstants, mascagnite::forwardRates, mascagnite::reverseRates, mascagnite::mobileSpeciesFlag );

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact
