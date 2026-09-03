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
#include "constitutive/ionicStrength/SpeciatedIonicStrength.hpp"
#include "constitutive/activity/Bdot.hpp"
#include "constitutive/activity/Identity.hpp"

namespace hpcReact
{

namespace geochemistry
{
// turn off uncrustify to allow for better readability of the parameters
// *****UNCRUSTIFY-OFF******

// ################################## Ultramafic rxn set ##################################
namespace ultramafics
{

constexpr CArrayWrapper<signed char, 21, 20> stoichMatrix =
{ //      OH-   CO2(aq)      CO3--      Mg2OH+++   Mg4(OH)++++   MgOH+   Mg2CO3++   MgCO3(aq)   MgHCO3+   Mg(H3SiO4)2      MgH2SiO4     MgH3SiO4+   H2SiO4--   H3SiO4-      H4(H2SiO4)----   H6(H2SiO4)--   H+      HCO3-      Mg++   SiO2(aq)
    {    -1,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         -1,       0,        0,        0   }, //  OH- + H+ = H2O
    {     0,        -1,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,          1,       1,        0,        0   }, //  CO2(aq) + H2O = HCO3- + H+
    {     0,         0,      -1,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         -1,       1,        0,        0   }, //  CO3-- + H+ = HCO3-
    {     0,         0,       0,         -1,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         -1,       0,        2,        0   }, //  Mg2OH+++ + H+ = 2Mg++ + H2O
    {     0,         0,       0,          0,         -1,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         -4,       0,        4,        0   }, //  Mg4(OH)++++ + 4H+ = 4Mg++ + 4H2O
    {     0,         0,       0,          0,          0,      -1,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         -1,       0,        1,        0   }, //  MgOH+ + H+ = Mg++ + H2O
    {     0,         0,       0,          0,          0,       0,         -1,           0,         0,             0,          0,           0,          0,         0,            0,              0,         -1,       1,        2,        0   }, //  Mg2CO3++ + H+ = 2Mg++ + HCO3-
    {     0,         0,       0,          0,          0,       0,          0,          -1,         0,             0,          0,           0,          0,         0,            0,              0,         -1,       1,        1,        0   }, //  MgCO3(aq) + H+ = Mg++ + HCO3-
    {     0,         0,       0,          0,          0,       0,          0,           0,        -1,             0,          0,           0,          0,         0,            0,              0,          0,       1,        1,        0   }, //  MgHCO3+ = Mg++ + HCO3-
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,            -1,          0,           0,          0,         0,            0,              0,         -2,       0,        1,        2   }, //  Mg(H3SiO4)2 + 2H+ = Mg++ + 2SiO2(aq) + 4H2O
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,         -1,           0,          0,         0,            0,              0,         -2,       0,        1,        1   }, //  MgH2SiO4 + 2H+ = Mg++ + SiO2(aq) + 2H2O
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,          -1,          0,         0,            0,              0,         -1,       0,        1,        1   }, //  MgH3SiO4+ + H+ = Mg++ + SiO2(aq) + 2H2O
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,         -1,         0,            0,              0,         -2,       0,        0,        1   }, //  H2SiO4-- + 2H+ = SiO2(aq) + 2H2O
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,        -1,            0,              0,         -1,       0,        0,        1   }, //  H3SiO4- + H+ = SiO2(aq) + 2H2O
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,           -1,              0,         -4,       0,        0,        4   }, //  H4(H2SiO4)---- + 4H+ = 4SiO2(aq) + 8H2O
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,             -1,         -2,       0,        0,        4   }, //  H6(H2SiO4)-- + 2H+ = 4SiO2 + 8H2O
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         -4,       0,        2,        1   }, //  Mg2SiO4(s) + 4H+ = 2Mg++ + SiO2(aq) + 2H2O (kinetic)
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         -1,       1,        1,        0   }, //  MgCO3(s) + H+ = Mg++ + HCO3- (kinetic)
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,          0,       0,        0,        1   }, //  SiO2(s) = SiO2(aq) (kinetic)
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         -6,       0,        3,        2   }, //  Mg3Si2O5(OH)4(s) + 6H+ = 3Mg++ + 2SiO2(aq) + 5H2O (kinetic)
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         -2,       0,        1,        0   }  //  Mg(OH)2(s) + 2H+ = Mg++ + 2H2O (kinetic)
  };

// 2Mg2SiO4 + 3H2O → Mg3Si2O5(OH)4 + Mg(OH)2 Serpentinization reaction  

constexpr CArrayWrapper<double, 21> equilibriumConstants = 
  { 
    9.89E+13,   //  OH- + H+ = H2O         
    4.42E-07,   //  CO2(aq) + H2O = HCO3- + H+  
    2.23E+10,   //  CO3-- + H+ = HCO3-       
    2.32E+13,   //  Mg2OH+++ + H+ = 2Mg++ + H2O
    4.47E+39,   //  Mg4(OH)++++ + 4H+ = 4Mg++ + 4H2O
    6.18E+11,   //  MgOH+ + H+ = Mg++ + H2O
    7.66E+06,   //  Mg2CO3++ + H+ = 2Mg++ + HCO3-
    2.67E+07,   //  MgCO3 + H+ = Mg++ + HCO3-
    9.77E-02,   //  MgHCO3+ = Mg++ + HCO3-
    3.45E+14,   //  Mg(H3SiO4)2 + 2H+ = Mg++ + 2SiO2(aq) + 4H2O
    9.49E+16,   //  MgH2SiO4 + 2H+ = Mg++ + SiO2(aq) + 2H2O
    1.96E+08,   //  MgH3SiO4+ + H+ = Mg++ + SiO2(aq) + 2H2O
    8.08E+22,   //  H2SiO4-- + 2H+ = SiO2(aq) + 2H2O
    6.44E+09,   //  H3SiO4- + H+ = SiO2(aq) + 2H2O
    5.39E+35,   //  H4(H2SiO4)---- + 4H+ = 4SiO2(aq) + 8H2O
    2.72E+13,   //  H6(H2SiO4)-- + 2H+ = 4SiO2 + 8H2O
    1.40E+28,   //  Mg2SiO4 + 4H+ = 2Mg++ + SiO2(aq) + 2H2O
    2.73E+02,   //  MgCO3 + H+ = Mg++ + HCO3-
    1.93E-03,   //  SiO2 = SiO2(aq)
    3.54E+31,   //  Mg3Si2O5(OH)4 + 6H+ = 3Mg++ + 2SiO2(aq) + 5H2O
    2.75E+16    //  Mg(OH)2 + 2H+ = Mg++ + 2H2O
  };

constexpr CArrayWrapper<double, 21> forwardRates = 
  { 
    1.00E+10,   //  OH- + H+ = H2O         
    1.00E+10,   //  CO2(aq) + H2O = HCO3- + H+  
    1.00E+10,   //  CO3-- + H+ = HCO3-       
    1.00E+10,   //  Mg2OH+++ + H+ = 2Mg++ + H2O
    1.00E+10,   //  Mg4(OH)++++ + 4H+ = 4Mg++ + 4H2O
    1.00E+10,   //  MgOH+ + H+ = Mg++ + H2O
    1.00E+10,   //  Mg2CO3++ + H+ = 2Mg++ + HCO3-
    1.00E+10,   //  MgCO3 + H+ = Mg++ + HCO3-
    1.00E+10,   //  MgHCO3+ = Mg++ + HCO3-
    1.00E+10,   //  Mg(H3SiO4)2 + 2H+ = Mg++ + 2SiO2(aq) + 4H2O
    1.00E+10,   //  MgH2SiO4 + 2H+ = Mg++ + SiO2(aq) + 2H2O
    1.00E+10,   //  MgH3SiO4+ + H+ = Mg++ + SiO2(aq) + 2H2O
    1.00E+10,   //  H2SiO4-- + 2H+ = SiO2(aq) + 2H2O
    1.00E+10,   //  H3SiO4- + H+ = SiO2(aq) + 2H2O
    1.00E+10,   //  H4(H2SiO4)---- + 4H+ = 4SiO2(aq) + 8H2O
    1.00E+10,   //  H6(H2SiO4)-- + 2H+ = 4SiO2 + 8H2O
    2.29E-11,   //  Mg2SiO4 + 4H+ = 2Mg++ + SiO2(aq) + 2H2O
    4.57E-10,   //  MgCO3 + H+ = Mg++ + HCO3-
    1.70E-13,   //  SiO2 = SiO2(aq)
    1.00E-12,   //  Mg3Si2O5(OH)4 + 6H+ = 3Mg++ + 2SiO2(aq) + 5H2O
    5.75E-09    //  Mg(OH)2 + 2H+ = Mg++ + 2H2O
  };

constexpr CArrayWrapper<double, 21> reverseRates = 
  { 
    1.00E+10,   //  OH- + H+ = H2O         
    1.00E+10,   //  CO2(aq) + H2O = HCO3- + H+  
    1.00E+10,   //  CO3-- + H+ = HCO3-       
    1.00E+10,   //  Mg2OH+++ + H+ = 2Mg++ + H2O
    1.00E+10,   //  Mg4(OH)++++ + 4H+ = 4Mg++ + 4H2O
    1.00E+10,   //  MgOH+ + H+ = Mg++ + H2O
    1.00E+10,   //  Mg2CO3++ + H+ = 2Mg++ + HCO3-
    1.00E+10,   //  MgCO3 + H+ = Mg++ + HCO3-
    1.00E+10,   //  MgHCO3+ = Mg++ + HCO3-
    1.00E+10,   //  Mg(H3SiO4)2 + 2H+ = Mg++ + 2SiO2(aq) + 4H2O
    1.00E+10,   //  MgH2SiO4 + 2H+ = Mg++ + SiO2(aq) + 2H2O
    1.00E+10,   //  MgH3SiO4+ + H+ = Mg++ + SiO2(aq) + 2H2O
    1.00E+10,   //  H2SiO4-- + 2H+ = SiO2(aq) + 2H2O
    1.00E+10,   //  H3SiO4- + H+ = SiO2(aq) + 2H2O
    1.00E+10,   //  H4(H2SiO4)---- + 4H+ = 4SiO2(aq) + 8H2O
    1.00E+10,   //  H6(H2SiO4)-- + 2H+ = 4SiO2 + 8H2O
    1.65E-39,   //  Mg2SiO4 + 4H+ = 2Mg++ + SiO2(aq) + 2H2O
    1.67E-12,   //  MgCO3 + H+ = Mg++ + HCO3-
    8.78E-11,   //  SiO2 = SiO2(aq)
    2.83E-44,   //  Mg3Si2O5(OH)4 + 6H+ = 3Mg++ + 2SiO2(aq) + 5H2O
    2.10E-25    //  Mg(OH)2 + 2H+ = Mg++ + 2H2O
  };

constexpr CArrayWrapper<int, 21> mobileSpeciesFlag = 
  { 
    1,   //  OH- + H+ = H2O         
    1,   //  CO2(aq) + H2O = HCO3- + H+  
    1,   //  CO3-- + H+ = HCO3-       
    1,   //  Mg2OH+++ + H+ = 2Mg++ + H2O
    1,   //  Mg4(OH)++++ + 4H+ = 4Mg++ + 4H2O
    1,   //  MgOH+ + H+ = Mg++ + H2O
    1,   //  Mg2CO3++ + H+ = 2Mg++ + HCO3-
    1,   //  MgCO3 + H+ = Mg++ + HCO3-
    1,   //  MgHCO3+ = Mg++ + HCO3-
    1,   //  Mg(H3SiO4)2 + 2H+ = Mg++ + 2SiO2(aq) + 4H2O
    1,   //  MgH2SiO4 + 2H+ = Mg++ + SiO2(aq) + 2H2O
    1,   //  MgH3SiO4+ + H+ = Mg++ + SiO2(aq) + 2H2O
    1,   //  H2SiO4-- + 2H+ = SiO2(aq) + 2H2O
    1,   //  H3SiO4- + H+ = SiO2(aq) + 2H2O
    1,   //  H4(H2SiO4)---- + 4H+ = 4SiO2(aq) + 8H2O
    1,   //  H6(H2SiO4)-- + 2H+ = 4SiO2 + 8H2O
    1,   //  Mg2SiO4 + 4H+ = 2Mg++ + SiO2(aq) + 2H2O
    1,   //  MgCO3 + H+ = Mg++ + HCO3-
    1,   //  SiO2 = SiO2(aq)
    1,   //  Mg3Si2O5(OH)4 + 6H+ = 3Mg++ + 2SiO2(aq) + 5H2O
    1    //  Mg(OH)2 + 2H+ = Mg++ + 2H2O
  };

// H2O coefficient, product-positive like the rows of stoichMatrix.
constexpr CArrayWrapper<signed char, 21> waterStoichiometry =
  {
     1,   //  OH- + H+ = H2O
    -1,   //  CO2(aq) + H2O = HCO3- + H+
     0,   //  CO3-- + H+ = HCO3-
     1,   //  Mg2OH+++ + H+ = 2Mg++ + H2O
     4,   //  Mg4(OH)++++ + 4H+ = 4Mg++ + 4H2O
     1,   //  MgOH+ + H+ = Mg++ + H2O
     0,   //  Mg2CO3++ + H+ = 2Mg++ + HCO3-
     0,   //  MgCO3(aq) + H+ = Mg++ + HCO3-
     0,   //  MgHCO3+ = Mg++ + HCO3-
     4,   //  Mg(H3SiO4)2 + 2H+ = Mg++ + 2SiO2(aq) + 4H2O
     2,   //  MgH2SiO4 + 2H+ = Mg++ + SiO2(aq) + 2H2O
     2,   //  MgH3SiO4+ + H+ = Mg++ + SiO2(aq) + 2H2O
     2,   //  H2SiO4-- + 2H+ = SiO2(aq) + 2H2O
     2,   //  H3SiO4- + H+ = SiO2(aq) + 2H2O
     8,   //  H4(H2SiO4)---- + 4H+ = 4SiO2(aq) + 8H2O
     8,   //  H6(H2SiO4)-- + 2H+ = 4SiO2 + 8H2O
     2,   //  Mg2SiO4(s) + 4H+ = 2Mg++ + SiO2(aq) + 2H2O
     0,   //  MgCO3(s) + H+ = Mg++ + HCO3-
     0,   //  SiO2(s) = SiO2(aq)
     5,   //  Mg3Si2O5(OH)4(s) + 6H+ = 3Mg++ + 2SiO2(aq) + 5H2O
     2    //  Mg(OH)2(s) + 2H+ = Mg++ + 2H2O
  };

// Activity model parameters.
//
// Charge z_i of each species, in the column order of stoichMatrix.
constexpr CArrayWrapper<double, 20> speciesCharge =
  {
    -1,   //  OH-
     0,   //  CO2(aq)
    -2,   //  CO3--
     3,   //  Mg2OH+++
     4,   //  Mg4(OH)++++
     1,   //  MgOH+
     2,   //  Mg2CO3++
     0,   //  MgCO3(aq)
     1,   //  MgHCO3+
     0,   //  Mg(H3SiO4)2
     0,   //  MgH2SiO4
     1,   //  MgH3SiO4+
    -2,   //  H2SiO4--
    -1,   //  H3SiO4-
    -4,   //  H4(H2SiO4)----
    -2,   //  H6(H2SiO4)--
     1,   //  H+
    -1,   //  HCO3-
     2,   //  Mg++
     0    //  SiO2(aq)
  };

// EQ3/6 B-dot parameters. Hard core diameters in ANGSTROM, the DHazero entries of data0.com.V8.R6,
// matching Mg4(OH)++++ -> Mg4(OH)4++++, H3SiO4- -> HSiO3-, and the polysilicates to their
// (H2SiO4)4 forms. MgOH+ is absent there and taken from data0.hmw. Mg2OH+++, Mg2CO3++ and
// MgH3SiO4+ are in no EQ3/6 database here and get the 4.0 default, which data0.com also gives every
// other charged complex in this system. The two neutral Mg silicates never use theirs.
constexpr CArrayWrapper<double, 20> ionSizeEQ36 =
  {
    3.5,   //  OH-
    3.0,   //  CO2(aq)
    4.5,   //  CO3--
    4.0,   //  Mg2OH+++       (default)
    5.5,   //  Mg4(OH)++++
    4.0,   //  MgOH+          (data0.hmw)
    4.0,   //  Mg2CO3++       (default)
    3.0,   //  MgCO3(aq)
    4.0,   //  MgHCO3+
    3.0,   //  Mg(H3SiO4)2    (neutral, unused)
    3.0,   //  MgH2SiO4       (neutral, unused)
    4.0,   //  MgH3SiO4+      (default)
    4.0,   //  H2SiO4--
    4.0,   //  H3SiO4-
    4.0,   //  H4(H2SiO4)----
    4.0,   //  H6(H2SiO4)--
    9.0,   //  H+
    4.0,   //  HCO3-
    8.0,   //  Mg++
    3.0    //  SiO2(aq)
  };

// The single b of data0.com.V8.R6 at 25 C. EQ3/6 applies it to charged species only.
constexpr double bdotEQ36_25C = 0.0410;

constexpr CArrayWrapper<double, 20> bdotParametersEQ36 =
  {
    bdotEQ36_25C,   //  OH-
    0.0,            //  CO2(aq)
    bdotEQ36_25C,   //  CO3--
    bdotEQ36_25C,   //  Mg2OH+++
    bdotEQ36_25C,   //  Mg4(OH)++++
    bdotEQ36_25C,   //  MgOH+
    bdotEQ36_25C,   //  Mg2CO3++
    0.0,            //  MgCO3(aq)
    bdotEQ36_25C,   //  MgHCO3+
    0.0,            //  Mg(H3SiO4)2
    0.0,            //  MgH2SiO4
    bdotEQ36_25C,   //  MgH3SiO4+
    bdotEQ36_25C,   //  H2SiO4--
    bdotEQ36_25C,   //  H3SiO4-
    bdotEQ36_25C,   //  H4(H2SiO4)----
    bdotEQ36_25C,   //  H6(H2SiO4)--
    bdotEQ36_25C,   //  H+
    bdotEQ36_25C,   //  HCO3-
    bdotEQ36_25C,   //  Mg++
    0.0             //  SiO2(aq)
  };

// EQ3/6 'neutral ion type' column. 0 is neutralSpeciesType::standard, -1 is
// neutralSpeciesType::drummond, which data0.com.V8.R6 gives to CO2(aq) alone.
constexpr CArrayWrapper<signed char, 20> neutralSpeciesTypeEQ36 =
  {
     0,   //  OH-
    -1,   //  CO2(aq)
     0,   //  CO3--
     0,   //  Mg2OH+++
     0,   //  Mg4(OH)++++
     0,   //  MgOH+
     0,   //  Mg2CO3++
     0,   //  MgCO3(aq)
     0,   //  MgHCO3+
     0,   //  Mg(H3SiO4)2
     0,   //  MgH2SiO4
     0,   //  MgH3SiO4+
     0,   //  H2SiO4--
     0,   //  H3SiO4-
     0,   //  H4(H2SiO4)----
     0,   //  H6(H2SiO4)--
     0,   //  H+
     0,   //  HCO3-
     0,   //  Mg++
     0    //  SiO2(aq)
  };
}

  // "all equilibrium" has no meaning once the minerals are gone: the five dissolution reactions
  // have no species left to be secondary, so 21 equilibrium reactions cannot be formed from 20
  // species.
  using ultramaficSystemAllKineticType     = reactionsSystems::MixedReactionsParameters< double, int, signed char, 20, 21, 0 >;
  using ultramaficSystemType               = reactionsSystems::MixedReactionsParameters< double, int, signed char, 20, 21, 16 >;

  // The species count of an activity model must match that of the system it is applied to, so it is
  // taken from the system type rather than repeated as a literal.
  using ultramaficIonicStrengthType    = SpeciatedIonicStrength< double, int, ultramaficSystemType::numSpecies() >;
  using ultramaficActivityType         = Bdot< double, int, ultramaficIonicStrengthType >;
  using ultramaficIdentityActivityType = Identity< double, int, ultramaficIonicStrengthType >;

  constexpr ultramaficSystemAllKineticType     ultramaficSystemAllKinetic( ultramafics::stoichMatrix, ultramafics::equilibriumConstants, ultramafics::forwardRates, ultramafics::reverseRates, ultramafics::mobileSpeciesFlag, reactionsSystems::ReactionRateLawOption::Affinity, ultramafics::waterStoichiometry );
  constexpr ultramaficSystemType               ultramaficSystem( ultramafics::stoichMatrix, ultramafics::equilibriumConstants, ultramafics::forwardRates, ultramafics::reverseRates, ultramafics::mobileSpeciesFlag, reactionsSystems::ReactionRateLawOption::Affinity, ultramafics::waterStoichiometry );

  constexpr ultramaficActivityType::Params ultramaficActivityParamsEQ36 =
  {
    {ultramafics::speciesCharge},
    ultramafics::ionSizeEQ36,
    ultramafics::bdotParametersEQ36,
    ultramafics::neutralSpeciesTypeEQ36,
    ultramafics::bdotEQ36_25C // the single b the water activity assumes all solutes share
  };

  constexpr ultramaficIdentityActivityType::Params ultramaficIdentityActivityParams = { { ultramafics::speciesCharge } };

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact
