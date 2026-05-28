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

// ################################## Ultramafic rxn set ##################################
namespace ultramafics
{

constexpr CArrayWrapper<signed char, 35, 32> stoichMatrix =
{ //      OH- CO2(aq) CO3-- CaCO3(aq) CaCl+ CaCl2(aq) CaHCO3+ CaHSiO3+ CaOH+ CaSO4(aq) HSO4- HSiO3- MgCO3(aq) MgCl+ MgHCO3+ MgHSiO3+ MgOH+ MgSO4(aq) NaCO3- NaCl(aq) NaHCO3(aq) NaHSiO3(aq) NaOH(aq) NaSO4-   H+  Cl- HCO3- Ca++ Mg++ Na+ SO4-- SiO2(aq)
    {    -1,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -1,   0,    0,   0,   0,   0,    0,    0 }, // OH- + H+ = H2O
    {     0,    -1,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,    1,   0,    1,   0,   0,   0,    0,    0 }, // CO2(aq) + H2O = H+ + HCO3-
    {     0,     0,    -1,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -1,   0,    1,   0,   0,   0,    0,    0 }, // CO3-- + H+ = HCO3-
    {     0,     0,     0,     -1,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -1,   0,    1,   1,   0,   0,    0,    0 }, // CaCO3(aq) + H+ = Ca++ + HCO3-
    {     0,     0,     0,      0,     -1,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,    0,   1,    0,   1,   0,   0,    0,    0 }, // CaCl+ = Ca++ + Cl-
    {     0,     0,     0,      0,      0,      -1,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,    0,   2,    0,   1,   0,   0,    0,    0 }, // CaCl2(aq) = Ca++ + 2Cl-
    {     0,     0,     0,      0,      0,       0,      -1,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,    0,   0,    1,   1,   0,   0,    0,    0 }, // CaHCO3+ = Ca++ + HCO3-
    {     0,     0,     0,      0,      0,       0,       0,       -1,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -1,   0,    0,   1,   0,   0,    0,    1 }, // CaHSiO3+ + H+ = Ca++ + SiO2(aq)
    {     0,     0,     0,      0,      0,       0,       0,        0,    -1,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -1,   0,    0,   1,   0,   0,    0,    0 }, // CaOH+ + H+ = Ca++ + H2O
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,      -1,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,    0,   0,    0,   1,   0,   0,    1,    0 }, // CaSO4(aq) = Ca++ + SO4--
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,    -1,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,    1,   0,    0,   0,   0,   0,    1,    0 }, // HSO4- = H+ + SO4--
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,     -1,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -1,   0,    0,   0,   0,   0,    0,    1 }, // HSiO3- + H+ = SiO2(aq)
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,      -1,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -1,   0,    1,   0,   1,   0,    0,    0 }, // MgCO3(aq) + H+ = Mg++ + HCO3-
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,   -1,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,    0,   1,    0,   0,   1,   0,    0,    0 }, // MgCl+ = Mg++ + Cl-
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,     -1,        0,    0,       0,      0,      0,        0,         0,      0,      0,    0,   0,    1,   0,   1,   0,    0,    0 }, // MgHCO3+ = Mg++ + HCO3-
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,       -1,    0,       0,      0,      0,        0,         0,      0,      0,   -1,   0,    0,   0,   1,   0,    0,    1 }, // MgHSiO3+ + H+ = Mg++ + SiO2(aq)
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,   -1,       0,      0,      0,        0,         0,      0,      0,   -1,   0,    0,   0,   1,   0,    0,    0 }, // MgOH+ + H+ = Mg++ + H2O
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,      -1,      0,      0,        0,         0,      0,      0,    0,   0,    0,   0,   1,   0,    1,    0 }, // MgSO4(aq) = Mg++ + SO4--
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,     -1,      0,        0,         0,      0,      0,   -1,   0,    1,   0,   0,   1,    0,    0 }, // NaCO3- + H+ = Na+ + HCO3-
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,     -1,        0,         0,      0,      0,    0,   1,    0,   0,   0,   1,    0,    0 }, // NaCl(aq) = Na+ + Cl-
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,       -1,         0,      0,      0,    0,   0,    1,   0,   0,   1,    0,    0 }, // NaHCO3(aq) = Na+ + HCO3-
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,        -1,      0,      0,   -1,   0,    0,   0,   0,   1,    0,    1 }, // NaHSiO3(aq) + H+ = Na+ + SiO2(aq)
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,     -1,      0,   -1,   0,    0,   0,   0,   1,    0,    0 }, // NaOH(aq) + H+ = Na+ + H2O
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,     -1,    0,   0,    0,   0,   0,   1,    1,    0 }, // NaSO4- = Na+ + SO4--
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -4,   0,    0,   0,   2,   0,    0,    1 }, // Forsterite + 4H+ = 2Mg++ + SiO2(aq)
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -2,   0,    0,   0,   1,   0,    0,    1 }, // Enstatite + 2H+ = Mg++ + SiO2(aq) + H2O
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -6,   0,    0,   0,   3,   0,    0,    2 }, // Chrysotile + 6H+ = 3Mg++ + 2SiO2(aq)
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -4,   0,    0,   1,   1,   0,    0,    2 }, // Diopside + 4H+ = Ca++ + Mg++ + 2SiO2(aq)
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -1,   0,    1,   1,   0,   0,    0,    0 }, // Calcite + H+ = Ca++ + HCO3-
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -2,   0,    2,   1,   1,   0,    0,    0 }, // Dolomite-ord + 2H+ = Ca++ + Mg++ + 2HCO3-
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -1,   0,    1,   0,   1,   0,    0,    0 }, // Magnesite + H+ = Mg++ + HCO3-
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -6,   0,    0,   0,   3,   0,    0,    4 }, // Talc + 6H+ = 3Mg++ + 4SiO2(aq)
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,  -14,   0,    0,   2,   5,   0,    0,    8 }, // Tremolite + 14H+ = 2Ca++ + 5Mg++ + 8SiO2(aq)
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,   -2,   0,    0,   0,   1,   0,    0,    0 }, // Brucite + 2H+ = Mg++ + 2H2O
    {     0,     0,     0,      0,      0,       0,       0,        0,     0,       0,     0,      0,       0,    0,      0,        0,    0,       0,      0,      0,        0,         0,      0,      0,    0,   0,    0,   0,   0,   0,    0,    1 }, // Quartz = SiO2(aq)
};

constexpr CArrayWrapper<double, 35> equilibriumConstants =
{
  1.06E+13, // OH- + H+ = H2O
  5.39E-07, // CO2(aq) + H2O = H+ + HCO3-
  1.35E+10, // CO3-- + H+ = HCO3-
  2.85E+06, // CaCO3(aq) + H+ = Ca++ + HCO3-
  1.24E+00, // CaCl+ = Ca++ + Cl-
  4.31E+00, // CaCl2(aq) = Ca++ + 2Cl-
  6.96E-02, // CaHCO3+ = Ca++ + HCO3-
  1.30E+08, // CaHSiO3+ + H+ = Ca++ + SiO2(aq)
  2.61E+11, // CaOH+ + H+ = Ca++ + H2O
  5.44E-03, // CaSO4(aq) = Ca++ + SO4--
  3.66E-03, // HSO4- = H+ + SO4--
  1.74E+09, // HSiO3- + H+ = SiO2(aq)
  8.50E+06, // MgCO3(aq) + H+ = Mg++ + HCO3-
  1.10E+00, // MgCl+ = Cl- + Mg++
  6.91E-02, // MgHCO3+ = HCO3- + Mg++
  8.78E+07, // MgHSiO3+ + H+ = Mg++ + SiO2(aq)
  3.17E+10, // MgOH+ + H+ = Mg++ + H2O
  4.05E-03, // MgSO4(aq) = Mg++ + SO4--
  8.27E+09, // NaCO3- + H+ = HCO3- + Na+
  4.48E+00, // NaCl(aq) = Cl- + Na+
  1.19E+00, // NaHCO3(aq) = HCO3- + Na+
  6.03E+07, // NaHSiO3(aq) + H+ = Na+ + SiO2(aq)
  1.62E+13, // NaOH(aq) + H+ = Na+ + H2O
  1.44E-01, // NaSO4- = Na+ + SO4--
  1.38E+24, // Forsterite + 4H+ = 2Mg++ + SiO2(aq)       (acid-neutral)
  7.20E+09, // Enstatite + 2H+ = Mg++ + SiO2(aq)         (acid-neutral)
  1.77E+27, // Chrysotile + 6H+ = 3Mg++ + 2SiO2(aq)      (neutral-base)
  4.02E+18, // Diopside + 4H+ = Ca++ + Mg++ + 2SiO2(aq)  (acid-neutral)
  2.15E+01, // Calcite + H+ = Ca++ + HCO3-               (acid-neutral-base)
  2.14E+01, // Dolomite-ord + 2H+ = Ca++ + Mg++ + 2HCO3- (acid-neutral-base)
  2.74E+01, // Magnesite + H+ = Mg++ + HCO3-             (acid-neutral-base)
  4.10E+18, // Talc + 6H+ = 3Mg++ + 4SiO2(aq)            (neutral)
  1.49E+54, // Tremolite + 14H+ = 2Ca++ + 5Mg++ + 8SiO2(aq) (acid-neutral)
  1.85E+14, // Brucite + 2H+ = Mg++ + 2H2O                  (acid-neutral)
  4.413E-04, // Quartz = SiO2(aq)                         (neutral-base)
};

constexpr CArrayWrapper<double, 35> forwardRates = 
  { 
  1.00E+10, // OH- + H+ = H2O
  1.00E+10, // CO2(aq) + H2O = H+ + HCO3-
  1.00E+10, // CO3-- + H+ = HCO3-
  1.00E+10, // CaCO3(aq) + H+ = Ca++ + HCO3-
  1.00E+10, // CaCl+ = Ca++ + Cl-
  1.00E+10, // CaCl2(aq) = Ca++ + 2Cl-
  1.00E+10, // CaHCO3+ = Ca++ + HCO3-
  1.00E+10, // CaHSiO3+ + H+ = Ca++ + SiO2(aq)
  1.00E+10, // CaOH+ + H+ = Ca++ + H2O
  1.00E+10, // CaSO4(aq) = Ca++ + SO4--
  1.00E+10, // HSO4- = H+ + SO4--
  1.00E+10, // HSiO3- + H+ = SiO2(aq)
  1.00E+10, // MgCO3(aq) + H+ = Mg++ + HCO3-
  1.00E+10, // MgCl+ = Cl- + Mg++
  1.00E+10, // MgHCO3+ = HCO3- + Mg++
  1.00E+10, // MgHSiO3+ + H+ = Mg++ + SiO2(aq)
  1.00E+10, // MgOH+ + H+ = Mg++ + H2O
  1.00E+10, // MgSO4(aq) = Mg++ + SO4--
  1.00E+10, // NaCO3- + H+ = HCO3- + Na+
  1.00E+10, // NaCl(aq) = Cl- + Na+
  1.00E+10, // NaHCO3(aq) = HCO3- + Na+
  1.00E+10, // NaHSiO3(aq) + H+ = Na+ + SiO2(aq)
  1.00E+10, // NaOH(aq) + H+ = Na+ + H2O
  1.00E+10, // NaSO4- = Na+ + SO4--
  4.77248E-09, // Forsterite + 4H+ = 2Mg++ + SiO2(aq)       (acid-neutral)
  4.24713E-11, // Enstatite + 2H+ = Mg++ + SiO2(aq)         (acid-neutral)
  1.43652E-10, // Chrysotile + 6H+ = 3Mg++ + 2SiO2(aq)      (neutral-base)
  1.20686E-10, // Diopside + 4H+ = Ca++ + Mg++ + 2SiO2(aq)  (acid-neutral)
  7.58146E-06, // Calcite + H+ = Ca++ + HCO3-               (acid-neutral-base)
  1.57459E-06, // Dolomite-ord + 2H+ = Ca++ + Mg++ + 2HCO3- (acid-neutral-base)
  2.23745E-09, // Magnesite + H+ = Mg++ + HCO3-             (acid-neutral-base)
  1.70903E-11, // Talc + 6H+ = 3Mg++ + 4SiO2(aq)            (neutral)
  1.48167E-08, // Tremolite + 14H+ = 2Ca++ + 5Mg++ + 8SiO2(aq) (acid-neutral)
  9.83445E-08,  // Brucite + 2H+ = Mg++ + 2H2O               (acid-neutral)
  1.85363E-11, // Quartz = SiO2(aq)                         (neutral-base)
  };

constexpr CArrayWrapper<double, 35> reverseRates = 
  { 
  1.00E+10, // OH- + H+ = H2O
  1.00E+10, // CO2(aq) + H2O = H+ + HCO3-
  1.00E+10, // CO3-- + H+ = HCO3-
  1.00E+10, // CaCO3(aq) + H+ = Ca++ + HCO3-
  1.00E+10, // CaCl+ = Ca++ + Cl-
  1.00E+10, // CaCl2(aq) = Ca++ + 2Cl-
  1.00E+10, // CaHCO3+ = Ca++ + HCO3-
  1.00E+10, // CaHSiO3+ + H+ = Ca++ + SiO2(aq)
  1.00E+10, // CaOH+ + H+ = Ca++ + H2O
  1.00E+10, // CaSO4(aq) = Ca++ + SO4--
  1.00E+10, // HSO4- = H+ + SO4--
  1.00E+10, // HSiO3- + H+ = SiO2(aq)
  1.00E+10, // MgCO3(aq) + H+ = Mg++ + HCO3-
  1.00E+10, // MgCl+ = Cl- + Mg++
  1.00E+10, // MgHCO3+ = HCO3- + Mg++
  1.00E+10, // MgHSiO3+ + H+ = Mg++ + SiO2(aq)
  1.00E+10, // MgOH+ + H+ = Mg++ + H2O
  1.00E+10, // MgSO4(aq) = Mg++ + SO4--
  1.00E+10, // NaCO3- + H+ = HCO3- + Na+
  1.00E+10, // NaCl(aq) = Cl- + Na+
  1.00E+10, // NaHCO3(aq) = HCO3- + Na+
  1.00E+10, // NaHSiO3(aq) + H+ = Na+ + SiO2(aq)
  1.00E+10, // NaOH(aq) + H+ = Na+ + H2O
  1.00E+10, // NaSO4- = Na+ + SO4--
  3.46E-33, // Forsterite + 4H+ = 2Mg++ + SiO2(aq)       (acid-neutral)
  5.90E-21, // Enstatite + 2H+ = Mg++ + SiO2(aq)         (acid-neutral)
  8.12E-38, // Chrysotile + 6H+ = 3Mg++ + 2SiO2(aq)      (neutral-base)
  3.00E-29, // Diopside + 4H+ = Ca++ + Mg++ + 2SiO2(aq)  (acid-neutral)
  3.53E-07, // Calcite + H+ = Ca++ + HCO3-               (acid-neutral-base)
  7.36E-08, // Dolomite-ord + 2H+ = Ca++ + Mg++ + 2HCO3- (acid-neutral-base)
  8.17E-11, // Magnesite + H+ = Mg++ + HCO3-             (acid-neutral-base)
  4.17E-30, // Talc + 6H+ = 3Mg++ + 4SiO2(aq)            (neutral)
  9.94E-63, // Tremolite + 14H+ = 2Ca++ + 5Mg++ + 8SiO2(aq) (acid-neutral)
  5.32E-22,  // Brucite + 2H+ = Mg++ + 2H2O               (acid-neutral)
  4.20E-08, // Quartz = SiO2(aq)                         (neutral-base)
  };

constexpr CArrayWrapper<int, 35> mobileSpeciesFlag =
  {
    1, // OH- + H+ = H2O
    1, // CO2(aq) + H2O = H+ + HCO3-
    1, // CO3-- + H+ = HCO3-
    1, // CaCO3(aq) + H+ = Ca++ + HCO3-
    1, // CaCl+ = Ca++ + Cl-
    1, // CaCl2(aq) = Ca++ + 2Cl-
    1, // CaHCO3+ = Ca++ + HCO3-
    1, // CaHSiO3+ + H+ = Ca++ + SiO2(aq) + H2O
    1, // CaOH+ + H+ = Ca++ + H2O
    1, // CaSO4(aq) = Ca++ + SO4--
    1, // HSO4- = H+ + SO4--
    1, // HSiO3- + H+ = SiO2(aq) + H2O
    1, // MgCO3(aq) + H+ = Mg++ + HCO3-
    1, // MgCl+ = Mg++ + Cl-
    1, // MgHCO3+ = Mg++ + HCO3-
    1, // MgHSiO3+ + H+ = Mg++ + SiO2(aq) + H2O
    1, // MgOH+ + H+ = Mg++ + H2O
    1, // MgSO4(aq) = Mg++ + SO4--
    1, // NaCO3- + H+ = Na+ + HCO3-
    1, // NaCl(aq) = Na+ + Cl-
    1, // NaHCO3(aq) = Na+ + HCO3-
    1, // NaHSiO3(aq) + H+ = Na+ + SiO2(aq) + H2O
    1, // NaOH(aq) + H+ = Na+ + H2O
    1, // NaSO4- = Na+ + SO4--
    1, // Forsterite + 4H+ = 2Mg++ + SiO2(aq) + 2H2O
    1, // Enstatite + 2H+ = Mg++ + SiO2(aq) + H2O
    1, // Chrysotile + 6H+ = 3Mg++ + 2SiO2(aq) + 5H2O
    1, // Diopside + 4H+ = Ca++ + Mg++ + 2SiO2(aq) + 2H2O
    1, // Calcite + H+ = Ca++ + HCO3-
    1, // Dolomite-ord + 2H+ = Ca++ + Mg++ + 2HCO3-
    1, // Magnesite + H+ = Mg++ + HCO3-
    1, // Talc + 6H+ = 3Mg++ + 4SiO2(aq) + 4H2O
    1, // Tremolite + 14H+ = 2Ca++ + 5Mg++ + 8SiO2(aq) + 8H2O
    1, // Brucite + 2H+ = Mg++ + 2H2O
    1, // Quartz = SiO2(aq)
  };
}

  //using ultramaficSystemAllKineticType     = reactionsSystems::MixedReactionsParameters< double, int, signed char, 32, 34, 0 >;
  //using ultramaficSystemAllEquilibriumType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 32, 34, 34 >;
  using ultramaficSystemType               = reactionsSystems::MixedReactionsParameters< double, int, signed char, 32, 35, 24 >;

  //constexpr ultramaficSystemAllKineticType     ultramaficSystemAllKinetic( ultramafics::stoichMatrix, ultramafics::equilibriumConstants, ultramafics::forwardRates, ultramafics::reverseRates, ultramafics::mobileSpeciesFlag );
  //constexpr ultramaficSystemAllEquilibriumType ultramaficSystemAllEquilibrium( ultramafics::stoichMatrix, ultramafics::equilibriumConstants, ultramafics::forwardRates, ultramafics::reverseRates, ultramafics::mobileSpeciesFlag );
  constexpr ultramaficSystemType               ultramaficSystem( ultramafics::stoichMatrix, ultramafics::equilibriumConstants, ultramafics::forwardRates, ultramafics::reverseRates, ultramafics::mobileSpeciesFlag );

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact
