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

// ################################## Primary serpentinization (Mg-olivine, no Fe) ##################################
// Net reaction: 2 Mg2SiO4 + 3 H2O → Mg3Si2O5(OH)4 + Mg(OH)2
//
// Primary species: H+, Mg++, SiO2(aq)
// Minerals (kinetic): Mg2SiO4 (forsterite), Mg3Si2O5(OH)4 (serpentine), Mg(OH)2 (brucite)
//
// All thermodynamic constants derived from 'llnl.tdat' (EQ3/6, via Geochemists' Workbench).
// Mineral forward rate constants from Palandri & Kharaka (2004) neutral-mechanism values at 25°C.

namespace serpentinization
{

constexpr CArrayWrapper<signed char, 3, 3> stoichMatrix =
{ //    H+     Mg++   SiO2(aq)
    {   -4,     2,       1   }, //  Mg2SiO4 + 4H+ = 2Mg++ + SiO2(aq) + 2H2O
    {   -6,     3,       2   }, //  Mg3Si2O5(OH)4 + 6H+ = 3Mg++ + 2SiO2(aq) + 5H2O
    {   -2,     1,       0   }  //  Mg(OH)2 + 2H+ = Mg++ + 2H2O
};

constexpr CArrayWrapper<double, 3> equilibriumConstants =
{
    1.40E+28,   //  Mg2SiO4 + 4H+ = 2Mg++ + SiO2(aq) + 2H2O
    3.54E+31,   //  Mg3Si2O5(OH)4 + 6H+ = 3Mg++ + 2SiO2(aq) + 5H2O
    2.75E+16    //  Mg(OH)2 + 2H+ = Mg++ + 2H2O
};

constexpr CArrayWrapper<double, 3> forwardRates =
{
    2.29E-11,   //  Mg2SiO4        [mol/(m2*s)]
    1.00E-12,   //  Mg3Si2O5(OH)4  [mol/(m2*s)]
    5.75E-09    //  Mg(OH)2        [mol/(m2*s)]
};

// kr = kf / K_eq
constexpr CArrayWrapper<double, 3> reverseRates =
{
    1.65E-39,   //  Mg2SiO4
    2.83E-44,   //  Mg3Si2O5(OH)4
    2.09E-25    //  Mg(OH)2
};

constexpr CArrayWrapper<int, 3> mobileSpeciesFlag =
{
    1,   //  Mg2SiO4
    1,   //  Mg3Si2O5(OH)4
    1    //  Mg(OH)2
};

} // namespace serpentinization

// 3 primary species, 3 kinetic mineral reactions, 0 equilibrium reactions
using serpentinizationSystemType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 3, 3, 0 >;

constexpr serpentinizationSystemType serpentinizationSystem( serpentinization::stoichMatrix, serpentinization::equilibriumConstants, serpentinization::forwardRates, serpentinization::reverseRates, serpentinization::mobileSpeciesFlag );

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact
