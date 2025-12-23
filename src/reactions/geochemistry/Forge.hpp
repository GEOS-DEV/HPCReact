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
// *****UNCRUSTIFY-OFF******

namespace forge
{

   // Stoichiometric matrix [13 reactions × 23 species]
   // Columns 0–12: secondary species (must be -1 on diagonal)
   // Columns 13–22: primary species
constexpr CArrayWrapper< signed char, 19, 26 > soichMatrix = 
{// CaCO₃  CaHCO₃⁺ CaSO₄ CaCl⁺ CaCl₂ MgHCO₃⁺ MgCO₃ MgCl⁺ CO₂(aq) HSO₄⁻ KHSO₄ HSiO₃⁻ NaHSilO₃ NaCl  KCl  KSO₄⁻ | H⁺  Ca²⁺  Mg²⁺   Na⁺    K⁺   Al³⁺   HCO₃⁻  SO₄²⁻  Cl⁻  SiO₂(aq)
  { -1,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,     -1,   1,    0,     0,    0,    0,     1,     0,     0,   0  }, // CaCO₃(aq) + H⁺ ⇌ Ca²⁺ + HCO₃⁻
  {  0,    -1,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,      0,   1,    0,     0,    0,    0,     1,     0,     0,   0  }, // CaHCO₃⁺ ⇌ Ca²⁺ + HCO₃⁻
  {  0,     0,    -1,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,      0,   1,    0,     0,    0,    0,     0,     0,     1,   0  }, // CaSO₄ ⇌ Ca²⁺ + SO₄²⁻
  {  0,     0,     0,   -1,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,      0,   1,    0,     0,    0,    0,     0,     0,     1,   0  }, // CaCl⁺ ⇌ Ca²⁺ + Cl⁻
  {  0,     0,     0,    0,   -1,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,      0,   1,    0,     0,    0,    0,     0,     0,     2,   0  }, // CaCl₂ ⇌ Ca²⁺ + 2Cl⁻
  {  0,     0,     0,    0,    0,    -1,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,      0,   0,    1,     0,    0,    0,     1,     0,     0,   0  }, // MgHCO₃⁺ ⇌ Mg²⁺ + HCO₃⁻
  {  0,     0,     0,    0,    0,     0,    -1,    0,     0,     0,     0,     0,     0,     0,   0,    0,     -1,   0,    1,     0,    0,    0,     1,     0,     0,   0  }, // MgCO₃(aq) + H⁺⇌ Mg²⁺ + HCO₃⁻
  {  0,     0,     0,    0,    0,     0,     0,   -1,     0,     0,     0,     0,     0,     0,   0,    0,      0,   0,    1,     0,    0,    0,     0,     0,     1,   0  }, // MgCl⁺ ⇌ Mg²⁺ + Cl⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,    -1,     0,     0,     0,     0,     0,   0,    0,      1,   0,    0,     0,    0,    0,     1,     0,     0,   0  }, // CO₂(aq) + H₂O ⇌ H⁺ + HCO₃⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,    -1,     0,     0,     0,     0,   0,    0,      1,   0,    0,     0,    0,    0,     0,     1,     0,   0  }, // HSO₄⁻ ⇌ H⁺ + SO₄²⁻ 
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,    -1,     0,     0,     0,   0,    0,      1,   0,    0,     0,    1,    0,     0,     1,     0,   0  }, // KHSO₄ ⇌ H⁺ + K⁺ + SO₄²⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,    -1,     0,     0,   0,    0,      1,   0,    0,     0,    0,    0,     0,     0,     0,   1  }, // HSiO₃⁻ ⇌ H⁺ + SiO₂(aq)
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,    -1,     0,   0,    0,      1,   0,    0,     1,    0,    0,     0,     0,     0,   1  }, // NaHSilO₃ ⇌ H⁺ + Na⁺ + SiO₂(aq)
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,    -1,   0,    0,      0,   0,    0,     1,    0,    0,     0,     0,     1,   0  }, // NaCl ⇌ Na⁺ + Cl⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,  -1,    0,      0,   0,    0,     0,    1,    0,     0,     0,     1,   0  }, // KCl ⇌ K⁺ + Cl⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,   -1,      0,   0,    0,     0,    1,    0,     0,     1,     0,   0  }, // KSO₄⁻ ⇌ K⁺ + SO₄²⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,     -2,   1,    1,     0,    0,    0,     2,     0,     0,   0  }, // Dolomite: CaMg(CO₃)₂(s) + 2H⁺ ⇌ Ca²⁺ + Mg²⁺ + 2HCO₃⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,     -4,   0,    0,     0,    1,    1,     0,     0,     0,   3  }, // Microcline: KAlSi₃O₈(s) + 4H⁺ ⇌ Al³⁺ + K⁺ + 3SiO₂(aq)
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,     -4,   0,    0,     1,    0,    1,     0,     1,     0,   3  }  // Albite: NaAlSi₃O₈(s) + 4H⁺ ⇌ Al³⁺ + Na⁺ + 3SiO₂(aq)
};

// Must convert these. They should not be the log.
constexpr CArrayWrapper< double, 19 > equilibriumConstants =
{
   5.9636,   // CaCO₃(aq) + H⁺ ⇌ Ca²⁺ + HCO₃⁻
  -1.4181,   // CaHCO₃⁺ ⇌ Ca²⁺ + HCO₃⁻
  -2.5111,   // CaSO₄ ⇌ Ca²⁺ + SO₄²⁻
   0.3811,   // CaCl⁺ ⇌ Ca²⁺ + Cl⁻
   0.3811,   // CaCl₂ ⇌ Ca²⁺ + 2Cl⁻ (approximate, same source)
  -1.4355,   // MgHCO₃⁺ ⇌ Mg²⁺ + HCO₃⁻
  -6.5632,   // MgCO₃(aq) + H⁺ ⇌ Mg²⁺ + HCO₃⁻
  -0.1820,   // MgCl⁺ ⇌ Mg²⁺ + Cl⁻
  -6.3882,   // CO₂(aq) + H₂O ⇌ H⁺ + HCO₃⁻
  -3.0020,   // HSO₄⁻ ⇌ H⁺ + SO₄²⁻
  -2.2935,   // KHSO₄ ⇌ H⁺ + K⁺ + SO₄²⁻
  -9.0844,   // HSiO₃⁻ ⇌ H⁺ + SiO₂(aq)
  -7.8291,   // NaHSilO₃ ⇌ H⁺ + Na⁺ + SiO₂(aq)
   0.4730,   // NaCl ⇌ Na⁺ + Cl⁻
   0.9240,   // KCl ⇌ K⁺ + Cl⁻
  -1.1946,   // KSO₄⁻ ⇌ K⁺ + SO₄²⁻
   0.0944,   // Dolomite: CaMg(CO₃)₂(s) + 2H⁺ ⇌ Ca²⁺ + Mg²⁺ + 2HCO₃⁻
  -1.8683,   // Microcline: KAlSi₃O₈(s) + 4H⁺ ⇌ Al³⁺ + K⁺ + 3SiO₂(aq)
   0.2236    // Albite: NaAlSi₃O₈(s) + 4H⁺ ⇌ Al³⁺ + Na⁺ + 3SiO₂(aq)
};

constexpr CArrayWrapper< double, 19 > fwRateConstant =
{
   0.0,   // CaCO₃(aq) + H⁺ ⇌ Ca²⁺ + HCO₃⁻
   0.0,   // CaHCO₃⁺ ⇌ Ca²⁺ + HCO₃⁻
   0.0,   // CaSO₄ ⇌ Ca²⁺ + SO₄²⁻
   0.0,   // CaCl⁺ ⇌ Ca²⁺ + Cl⁻
   0.0,   // CaCl₂ ⇌ Ca²⁺ + 2Cl⁻ (approximate, same source)
   0.0,   // MgHCO₃⁺ ⇌ Mg²⁺ + HCO₃⁻
   0.0,   // MgCO₃(aq) + H⁺ ⇌ Mg²⁺ + HCO₃⁻
   0.0,   // MgCl⁺ ⇌ Mg²⁺ + Cl⁻
   0.0,   // CO₂(aq) + H₂O ⇌ H⁺ + HCO₃⁻
   0.0,   // HSO₄⁻ ⇌ H⁺ + SO₄²⁻
   0.0,   // KHSO₄ ⇌ H⁺ + K⁺ + SO₄²⁻
   0.0,   // HSiO₃⁻ ⇌ H⁺ + SiO₂(aq)
   0.0,   // NaHSilO₃ ⇌ H⁺ + Na⁺ + SiO₂(aq)
   0.0,   // NaCl ⇌ Na⁺ + Cl⁻
   0.0,   // KCl ⇌ K⁺ + Cl⁻
   0.0,   // KSO₄⁻ ⇌ K⁺ + SO₄²⁻
   1.0,   // Dolomite: CaMg(CO₃)₂(s) + 2H⁺ ⇌ Ca²⁺ + Mg²⁺ + 2HCO₃⁻
   1.0,   // Microcline: KAlSi₃O₈(s) + 4H⁺ ⇌ Al³⁺ + K⁺ + 3SiO₂(aq)
   1.0    // Albite: NaAlSi₃O₈(s) + 4H⁺ ⇌ Al³⁺ + Na⁺ + 3SiO₂(aq)
};


constexpr CArrayWrapper< double, 19 > reverseRateConstant =
{
   0.0,   // CaCO₃(aq) + H⁺ ⇌ Ca²⁺ + HCO₃⁻
   0.0,   // CaHCO₃⁺ ⇌ Ca²⁺ + HCO₃⁻
   0.0,   // CaSO₄ ⇌ Ca²⁺ + SO₄²⁻
   0.0,   // CaCl⁺ ⇌ Ca²⁺ + Cl⁻
   0.0,   // CaCl₂ ⇌ Ca²⁺ + 2Cl⁻ (approximate, same source)
   0.0,   // MgHCO₃⁺ ⇌ Mg²⁺ + HCO₃⁻
   0.0,   // MgCO₃(aq) + H⁺ ⇌ Mg²⁺ + HCO₃⁻
   0.0,   // MgCl⁺ ⇌ Mg²⁺ + Cl⁻
   0.0,   // CO₂(aq) + H₂O ⇌ H⁺ + HCO₃⁻
   0.0,   // HSO₄⁻ ⇌ H⁺ + SO₄²⁻
   0.0,   // KHSO₄ ⇌ H⁺ + K⁺ + SO₄²⁻
   0.0,   // HSiO₃⁻ ⇌ H⁺ + SiO₂(aq)
   0.0,   // NaHSilO₃ ⇌ H⁺ + Na⁺ + SiO₂(aq)
   0.0,   // NaCl ⇌ Na⁺ + Cl⁻
   0.0,   // KCl ⇌ K⁺ + Cl⁻
   0.0,   // KSO₄⁻ ⇌ K⁺ + SO₄²⁻
   1.0,   // Dolomite: CaMg(CO₃)₂(s) + 2H⁺ ⇌ Ca²⁺ + Mg²⁺ + 2HCO₃⁻
   1.0,   // Microcline: KAlSi₃O₈(s) + 4H⁺ ⇌ Al³⁺ + K⁺ + 3SiO₂(aq)
   1.0    // Albite: NaAlSi₃O₈(s) + 4H⁺ ⇌ Al³⁺ + Na⁺ + 3SiO₂(aq)
};

constexpr CArrayWrapper< int, 19 > mobileSpeciesFlag =
{
   1,   // CaCO₃(aq) + H⁺ ⇌ Ca²⁺ + HCO₃⁻
   1,   // CaHCO₃⁺ ⇌ Ca²⁺ + HCO₃⁻
   1,   // CaSO₄ ⇌ Ca²⁺ + SO₄²⁻
   1,   // CaCl⁺ ⇌ Ca²⁺ + Cl⁻
   1,   // CaCl₂ ⇌ Ca²⁺ + 2Cl⁻ (approximate, same source)
   1,   // MgHCO₃⁺ ⇌ Mg²⁺ + HCO₃⁻
   1,   // MgCO₃(aq) + H⁺ ⇌ Mg²⁺ + HCO₃⁻
   1,   // MgCl⁺ ⇌ Mg²⁺ + Cl⁻
   1,   // CO₂(aq) + H₂O ⇌ H⁺ + HCO₃⁻
   1,   // HSO₄⁻ ⇌ H⁺ + SO₄²⁻
   1,   // KHSO₄ ⇌ H⁺ + K⁺ + SO₄²⁻
   1,   // HSiO₃⁻ ⇌ H⁺ + SiO₂(aq)
   1,   // NaHSilO₃ ⇌ H⁺ + Na⁺ + SiO₂(aq)
   1,   // NaCl ⇌ Na⁺ + Cl⁻
   1,   // KCl ⇌ K⁺ + Cl⁻
   1,   // KSO₄⁻ ⇌ K⁺ + SO₄²⁻
   1,   // Dolomite: CaMg(CO₃)₂(s) + 2H⁺ ⇌ Ca²⁺ + Mg²⁺ + 2HCO₃⁻
   1,   // Microcline: KAlSi₃O₈(s) + 4H⁺ ⇌ Al³⁺ + K⁺ + 3SiO₂(aq)
   1    // Albite: NaAlSi₃O₈(s) + 4H⁺ ⇌ Al³⁺ + Na⁺ + 3SiO₂(aq)
};

}

using forgeSystemType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 26, 19, 16 >;


constexpr forgeSystemType forgeSystem( forge::soichMatrix, forge::equilibriumConstants, forge::fwRateConstant, forge::reverseRateConstant, forge::mobileSpeciesFlag );


// *****UNCRUSTIFY-ON******
}
}
