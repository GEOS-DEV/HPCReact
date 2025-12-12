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

   // Stoichiometric matrix [24 reactions × 29 species]
   // Columns 1–19: secondary species (must be -1 on diagonal)
   // Columns 20–29: primary species
constexpr CArrayWrapper< double, 24, 29 > stoichMatrix = 
{// CaCO₃  CaHCO₃⁺ CaSO₄ CaCl⁺ CaCl₂ MgHCO₃⁺ MgCO₃ MgCl⁺ CO₂(aq) HSO₄⁻ KHSO₄ HSiO₃⁻ NaHSilO₃ NaCl  KCl  KSO₄⁻ AlOH²⁺ Al(OH)₂ OH⁻  | H⁺  Ca²⁺  Mg²⁺   Na⁺    K⁺   Al³⁺   HCO₃⁻  SO₄²⁻  Cl⁻  SiO₂(aq)
  { -1,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,    -1,   1,    0,     0,    0,    0,     1,     0,     0,   0  }, // CaCO₃(aq) + H⁺ ⇌ Ca²⁺ + HCO₃⁻
  {  0,    -1,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,     0,   1,    0,     0,    0,    0,     1,     0,     0,   0  }, // CaHCO₃⁺ ⇌ Ca²⁺ + HCO₃⁻
  {  0,     0,    -1,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,     0,   1,    0,     0,    0,    0,     0,     0,     1,   0  }, // CaSO₄(aq) ⇌ Ca²⁺ + SO₄²⁻
  {  0,     0,     0,   -1,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,     0,   1,    0,     0,    0,    0,     0,     0,     1,   0  }, // CaCl⁺ ⇌ Ca²⁺ + Cl⁻
  {  0,     0,     0,    0,   -1,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,     0,   1,    0,     0,    0,    0,     0,     0,     2,   0  }, // CaCl₂(aq) ⇌ Ca²⁺ + 2Cl⁻
  {  0,     0,     0,    0,    0,    -1,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,     0,   0,    1,     0,    0,    0,     1,     0,     0,   0  }, // MgHCO₃⁺ ⇌ Mg²⁺ + HCO₃⁻
  {  0,     0,     0,    0,    0,     0,    -1,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,    -1,   0,    1,     0,    0,    0,     1,     0,     0,   0  }, // MgCO₃(aq) + H⁺⇌ Mg²⁺ + HCO₃⁻
  {  0,     0,     0,    0,    0,     0,     0,   -1,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,     0,   0,    1,     0,    0,    0,     0,     0,     1,   0  }, // MgCl⁺ ⇌ Mg²⁺ + Cl⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,    -1,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,     1,   0,    0,     0,    0,    0,     1,     0,     0,   0  }, // CO₂(aq) + H₂O ⇌ H⁺ + HCO₃⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,    -1,     0,     0,     0,     0,   0,    0,    0,     0,      0,     1,   0,    0,     0,    0,    0,     0,     1,     0,   0  }, // HSO₄⁻ ⇌ H⁺ + SO₄²⁻ 
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,    -1,     0,     0,     0,   0,    0,    0,     0,      0,     1,   0,    0,     0,    1,    0,     0,     1,     0,   0  }, // KHSO₄(aq) ⇌ H⁺ + K⁺ + SO₄²⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,    -1,     0,     0,   0,    0,    0,     0,      0,     1,   0,    0,     0,    0,    0,     0,     0,     0,   1  }, // HSiO₃⁻ ⇌ H⁺ + SiO₂(aq)
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,    -1,     0,   0,    0,    0,     0,      0,     1,   0,    0,     1,    0,    0,     0,     0,     0,   1  }, // NaHSiO₃(aq) ⇌ H⁺ + Na⁺ + SiO₂(aq)
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,    -1,   0,    0,    0,     0,      0,     0,   0,    0,     1,    0,    0,     0,     0,     1,   0  }, // NaCl(aq) ⇌ Na⁺ + Cl⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,  -1,    0,    0,     0,      0,     0,   0,    0,     0,    1,    0,     0,     0,     1,   0  }, // KCl(aq) ⇌ K⁺ + Cl⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,   -1,    0,     0,      0,     0,   0,    0,     0,    1,    0,     0,     1,     0,   0  }, // KSO₄⁻ ⇌ K⁺ + SO₄²⁻
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,   -1,     0,      0,    -1,   0,    0,     0,    0,    1,     0,     0,     0,   0  }, // AlOH²⁺ + H⁺ ⇌ Al³⁺ + H₂O
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,    -1,      0,    -2,   0,    0,    0,     0,    1,     0,     0,     0,   0  }, // Al(OH)₂⁺ + 2H⁺ ⇌ Al³⁺ + 2H₂O
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,     -1,    -1,   0,    0,     0,    0,    0,     0,     0,     0,   0  }, // OH⁻ + H⁺ ⇌ H₂O
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,     0,   0,    0,     0,    0,    0,     0,     0,     0,   1  }, // SiO2(s) ⇌ SiO2(aq)
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,   -10,   0,    3,     0,    1,    1,     0,     0,     0,   3  }, // KAlMg3Si3O10(OH)2(s) + 10H+ ⇌ Al3+ + K+ + 3Mg2+ + 3SiO2(aq) + 6H2O
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,    -8,   1,    0,     0,    0,    2,     0,     0,     0,   2  }, // CaAl2(SiO4)2(s) + 8H+ ⇌ Ca2+ + 2 Al3+ + 2 SiO2(aq) + 4H2O
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,    -4,   0,    0,     0,    1,    1,     0,     0,     0,   3  }, // KAlSi3O8(s) + 4H+ ⇌ Al3+ + K+ + 3 SiO2 + 2H2O
  {  0,     0,     0,    0,    0,     0,     0,    0,     0,     0,     0,     0,     0,     0,   0,    0,    0,     0,      0,    -6,   0,    0,     0,    0,    2,     0,     0,     0,   2  }  // Al2Si2O5(OH)4(s) + 6H+ ⇌ 2Al3+ + 2 SiO2 + 5H2O
};

// Must convert these. They should not be the log.
constexpr CArrayWrapper< double, 24 > equilibriumConstants =
{
   104351.8133, // CaCO₃(aq) + H⁺ ⇌ Ca²⁺ + HCO₃⁻
   3.9739, //3.98107E-03, // CaHCO₃⁺ ⇌ Ca²⁺ + HCO₃⁻
   0.3685, //3.6915E-04,  // CaSO₄ ⇌ Ca²⁺ + SO₄²⁻
   292.4726, //0.2930     // CaCl⁺ ⇌ Ca²⁺ + Cl⁻
   1.2236e+05, //0.1228     // CaCl₂ ⇌ Ca²⁺ + 2Cl⁻ (approximate, same source)
   3.7932, //0.0038     // MgHCO₃⁺ ⇌ Mg²⁺ + HCO₃⁻
   7.4559E+05,  // MgCO₃(aq) + H⁺ ⇌ Mg²⁺ + HCO₃⁻
   72.4693,  //0.0726    // MgCl⁺ ⇌ Mg²⁺ + Cl⁻
   6.3434e-05, //6.3548E-08  // CO₂(aq) + H₂O ⇌ H⁺ + HCO₃⁻
   0.0340, //3.4017E-05,  // HSO₄⁻ ⇌ H⁺ + SO₄²⁻
   36.9347, //3.7068E-05,  // KHSO₄ ⇌ H⁺ + K⁺ + SO₄²⁻
   6.8964e+11, // 6.9088E+08,  // HSiO₃⁻ ⇌ H⁺ + SiO₂(aq)
   4.5356e+13, //4.5520E+07,  // NaHSiO₃ ⇌ H⁺ + Na⁺ + SiO₂(aq)
   805.2479, //0.8067,      // NaCl ⇌ Na⁺ + Cl⁻
   1.6368e+03, //1.6398,      // KCl ⇌ K⁺ + Cl⁻
   12.0782, //0.0121,      // KSO₄⁻ ⇌ K⁺ + SO₄²⁻
   0.0147, //14.6994,     // AlOH²⁺ + H⁺ ⇌ Al³⁺ + H₂O
   0.0039, //3.8530E+03,  // Al(OH)₂⁺ + 2H⁺ ⇌ Al³⁺ + 2H₂O
   1.9282e+05, //1.9213E+11,  // OH⁻ + H⁺ ⇌ H₂O
   0.0036,      // SiO2(s) ⇌ SiO2(aq)
   6.6130E+15,  // KAlMg3Si3O10(OH)2(s) + 10H+ ⇌ Al3+ + K+ + 3Mg2+ + 3SiO2(aq) + 6H2O
   1.3047E+05,  // CaAl2(SiO4)2(s) + 8H+ ⇌ Ca2+ + 2 Al3+ + 2 SiO2(aq) + 4H2O
   1.7669E-04,  // KAlSi3O8(s) + 4H+ ⇌ Al3+ + K+ + 3 SiO2 + 2H2O
   3.1463E-05   // Al2Si2O5(OH)4(s) + 6H+ ⇌ 2Al3+ + 2 SiO2 + 5H2O
};

constexpr CArrayWrapper< double, 24 > fwRateConstant =
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
   0.0,   // AlOH²⁺ + H⁺ ⇌ Al³⁺ + H₂O
   0.0,   // Al(OH)₂⁺ + 2H⁺ ⇌ Al³⁺ + 2H₂O
   0.0,   // OH⁻ + H⁺ ⇌ H₂O
   2.7043E-08,  // SiO2(s) ⇌ SiO2(aq)
   3.0145E-11,  // KAlMg3Si3O10(OH)2(s) + 10H+ ⇌ Al3+ + K+ + 3Mg2+ + 3SiO2(aq) + 6H2O
   1.0801E-08,  // CaAl2(SiO4)2(s) + 8H+ ⇌ Ca2+ + 2 Al3+ + 2 SiO2(aq) + 4H2O
   1.1283e-10,  // KAlSi3O8(s) + 4H+ ⇌ Al3+ + K+ + 3 SiO2 + 2H2O
   1.8137e-12   // Al2Si2O5(OH)4(s) + 6H+ ⇌ 2Al3+ + 2 SiO2 + 5H2O
};


constexpr CArrayWrapper< double, 24 > reverseRateConstant =
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
   0.0,   // AlOH²⁺ + H⁺ ⇌ Al³⁺ + H₂O
   0.0,   // Al(OH)₂⁺ + 2H⁺ ⇌ Al³⁺ + 2H₂O
   0.0,   // OH⁻ + H⁺ ⇌ H₂O
   7.4240e-09, //7.5119E-06,  // SiO2(s) ⇌ SiO2(aq)
   4.5420e-21, //4.5584E-27,  // KAlMg3Si3O10(OH)2(s) + 10H+ ⇌ Al3+ + K+ + 3Mg2+ + 3SiO2(aq) + 6H2O
   8.2341e-05, //8.2785E-14,  // CaAl2(SiO4)2(s) + 8H+ ⇌ Ca2+ + 2 Al3+ + 2 SiO2(aq) + 4H2O
   6.3975e-10, //6.3858E-07,  // KAlSi3O8(s) + 4H+ ⇌ Al3+ + K+ + 3 SiO2 + 2H2O
   0.0574     //5.7645E-08   // Al2Si2O5(OH)4(s) + 6H+ ⇌ 2Al3+ + 2 SiO2 + 5H2O
};

constexpr CArrayWrapper< int, 24 > mobileSpeciesFlag =
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
   1,   // AlOH²⁺ + H⁺ ⇌ Al³⁺ + H₂O
   1,   // Al(OH)₂⁺ + 2H⁺ ⇌ Al³⁺ + 2H₂O
   1,   // OH⁻ + H⁺ ⇌ H₂O
   1,  // SiO2(s) ⇌ SiO2(aq)
   1,  // KAlMg3Si3O10(OH)2(s) + 10H+ ⇌ Al3+ + K+ + 3Mg2+ + 3SiO2(aq) + 6H2O
   1,  // CaAl2(SiO4)2(s) + 8H+ ⇌ Ca2+ + 2 Al3+ + 2 SiO2(aq) + 4H2O
   1,  // KAlSi3O8(s) + 4H+ ⇌ Al3+ + K+ + 3 SiO2 + 2H2O
   1   // Al2Si2O5(OH)4(s) + 6H+ ⇌ Al3+ + 2 SiO2 + 5H2O
};

}

using forgeSystemType = reactionsSystems::MixedReactionsParameters< double, int, int, 29, 24, 19 >;


constexpr forgeSystemType forgeSystem( forge::stoichMatrix, forge::equilibriumConstants, forge::fwRateConstant, forge::reverseRateConstant, forge::mobileSpeciesFlag );


// *****UNCRUSTIFY-ON******
}
}