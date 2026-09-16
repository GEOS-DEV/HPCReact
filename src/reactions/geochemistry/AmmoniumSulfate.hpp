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

// ################################## Minimal ammonium sulfate system ##################################
// Single kinetic reaction: (NH4)2SO4(s) = 2 NH4+ + SO4--   (mascagnite, positive rate = dissolution)
//
// Primary species: NH4+, SO4--
// Mineral (kinetic): (NH4)2SO4
//
// Aqueous speciation (NH3/NH4+, HSO4-/SO4--) is neglected: at pH 5-7 NH3 and HSO4- are below 1% of
// their parent ions. Activities are ideal (gamma = 1), so K is a conditional solubility product
// calibrated to the measured solubility s = 5.78 mol/kg at 25 C (43.3 wt%): K = 4 s^3 = 772.402.
// The rate constant is of the order of halite (Palandri & Kharaka, 2004) and is a calibration knob.
//
// Equilibrium: Q/K = [NH4+]^2 [SO4--] / 772.402
// Saturated 2:1 solution: [SO4--] = 5.78 mol/kg, [NH4+] = 11.56 mol/kg gives Q/K = 1 exactly

namespace ammoniumSulfate
{

constexpr CArrayWrapper<signed char, 1, 2> stoichMatrix =
{ //   NH4+   SO4--
    {    2,     1   }  //  (NH4)2SO4(s) = 2 NH4+ + SO4--
};

constexpr CArrayWrapper<double, 1> equilibriumConstants =
{
    7.724022E+02   //  (NH4)2SO4(s) = 2 NH4+ + SO4--
};

constexpr CArrayWrapper<double, 1> forwardRates =
{
    1.0E-01   //  (NH4)2SO4  [mol/(m2*s)]
};

// kr = kf / K_eq
constexpr CArrayWrapper<double, 1> reverseRates =
{
    1.0E-01 / 7.724022E+02   //  (NH4)2SO4
};

constexpr CArrayWrapper<int, 1> mobileSpeciesFlag =
{
    1          //  (NH4)2SO4
};

} // namespace ammoniumSulfate

// 2 primary species, 1 kinetic mineral reaction, 0 equilibrium reactions
using ammoniumSulfateSystemType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 2, 1, 0 >;

constexpr ammoniumSulfateSystemType ammoniumSulfateSystem( ammoniumSulfate::stoichMatrix, ammoniumSulfate::equilibriumConstants, ammoniumSulfate::forwardRates, ammoniumSulfate::reverseRates, ammoniumSulfate::mobileSpeciesFlag );

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact
