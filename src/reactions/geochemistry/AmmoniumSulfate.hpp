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
//
// Temperature dependence. AS dissolution is endothermic, so solubility rises steeply with T and a solution
// made up warm supersaturates on cooling. K(T) follows
//     ln K(T) = ln K(T_ref) + B (1/T - 1/T_ref) + C ln(T/T_ref),   T_ref = 298.15 K
// least-squares fitted to measured solubility over 0-100 C with K(T_ref) pinned at 4 s^3 = 772.402, so the
// 25 C calibration point is exact and a fluid saturated at 25 C is still exactly Q/K = 1. Worst error over
// the range is 1.7%. A single-enthalpy van't Hoff form was rejected: accurate to 25 C but -33% at 100 C.
// Cooling a fluid saturated at T_prep down to 25 C gives Q/K = 1.18 (40 C), 1.50 (60 C), 1.94 (80 C),
// 2.53 (100 C).

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

// ln K(T) coefficients, see the note above. Zero would mean a temperature-independent K.
constexpr CArrayWrapper<double, 1> lnEqConstCoeffB =
{
    2.682837E+03   //  (NH4)2SO4  [K]
};

constexpr CArrayWrapper<double, 1> lnEqConstCoeffC =
{
    1.219289E+01   //  (NH4)2SO4  [-]
};

constexpr double referenceTemperature = 298.15;   // [K]

} // namespace ammoniumSulfate

// 2 primary species, 1 kinetic mineral reaction, 0 equilibrium reactions
using ammoniumSulfateSystemType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 2, 1, 0 >;

constexpr ammoniumSulfateSystemType ammoniumSulfateSystem( ammoniumSulfate::stoichMatrix, ammoniumSulfate::equilibriumConstants, ammoniumSulfate::forwardRates, ammoniumSulfate::reverseRates, ammoniumSulfate::mobileSpeciesFlag,
                                                          1, hpcReact::constants::waterDensity,
                                                          ammoniumSulfate::lnEqConstCoeffB, ammoniumSulfate::lnEqConstCoeffC, ammoniumSulfate::referenceTemperature );

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact
