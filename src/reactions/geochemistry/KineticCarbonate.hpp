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

// ################################## Minimal kinetic carbonate system ##################################
// Single kinetic reaction: CaCO3(s) + H+ = Ca+2 + HCO3-
//
// Primary species: H+, Ca+2, HCO3-
// Mineral (kinetic): CaCO3 (calcite)
//
// Thermodynamic constant from 'llnl.tdat' (EQ3/6, via Geochemists' Workbench).
// Mineral rate constant from Plummer et al. (1978) / Palandri & Kharaka (2004) neutral mechanism at 25°C.
//
// Equilibrium: Q/K = [Ca+2][HCO3-]/[H+] = 51.6
// At pH=8: [Ca+2] = [HCO3-] = sqrt(51.6 * 1e-8) = 7.18e-4 M (zero-rate initial condition)

namespace kineticCarbonate
{

constexpr CArrayWrapper<signed char, 1, 3> stoichMatrix =
{ //    H+     Ca+2   HCO3-
    {   -1,     1,       1   }  //  CaCO3(s) + H+ = Ca+2 + HCO3-
};

constexpr CArrayWrapper<double, 1> equilibriumConstants =
{
    4.9E+01   //  CaCO3(s) + H+ = Ca+2 + HCO3-
};

constexpr CArrayWrapper<double, 1> forwardRates =
{
    1.47E-06   //  CaCO3  [mol/(m2*s)]
};

// kr = kf / K_eq
constexpr CArrayWrapper<double, 1> reverseRates =
{
    3.00E-08   //  CaCO3
};

constexpr CArrayWrapper<int, 1> mobileSpeciesFlag =
{
    1          //  CaCO3
};

} // namespace kineticCarbonate

// 3 primary species, 1 kinetic mineral reaction, 0 equilibrium reactions
using kineticCarbonateSystemType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 3, 1, 0 >;

constexpr kineticCarbonateSystemType kineticCarbonateSystem( kineticCarbonate::stoichMatrix, kineticCarbonate::equilibriumConstants, kineticCarbonate::forwardRates, kineticCarbonate::reverseRates, kineticCarbonate::mobileSpeciesFlag );

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact