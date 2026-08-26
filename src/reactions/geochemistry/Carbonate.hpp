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

namespace carbonate
{

constexpr CArrayWrapper<signed char, 10, 17> stoichMatrix = 
  { //   OH-    CO2  CO3-2  CaHCO3+   CaSO4  CaCl+  CaCl2  MgSO4   NaSO4- CaCO3  H+  HCO3-  Ca+2    SO4-2    Cl-    Mg+2  Na+
    {    -1,     0,     0,      0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0  }, //     OH- + H+ = H2O         
    {     0,    -1,     0,      0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0,     0,     0  }, //    CO2 + H2O = H+ + HCO3-  
    {     0,     0,    -1,      0,     0,     0,     0,     0,     0,     0,    -1,     1,     0,     0,     0,     0,     0  }, //   CO3-2 + H+ = HCO3-       
    {     0,     0,     0,     -1,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0,     0  }, //      CaHCO3+ = Ca+2 + HCO3-
    {     0,     0,     0,      0,    -1,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0  }, //        CaSO4 = Ca+2 + SO4-2
    {     0,     0,     0,      0,     0,    -1,     0,     0,     0,     0,     0,     0,     1,     0,     1,     0,     0  }, //        CaCl+ = Ca+2 + Cl-  
    {     0,     0,     0,      0,     0,     0,    -1,     0,     0,     0,     0,     0,     1,     0,     2,     0,     0  }, //        CaCl2 = Ca+2 + 2Cl- 
    {     0,     0,     0,      0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     1,     0,     1,     0  }, //        MgSO4 = Mg+2 + SO4-2
    {     0,     0,     0,      0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     1,     0,     0,     1  }, //       NaSO4- = Na+ + SO4-2
    {     0,     0,     0,      0,     0,     0,     0,     0,     0,    -1,    -1,     1,     1,     0,     0,     0,     0  }  //   CaCO3(s) + H+ = Ca+2 + HCO3- (kinetic)
  };

constexpr CArrayWrapper<signed char, 10, 16> stoichMatrixNosolid = 
  { //   OH-    CO2  CO3-2  CaHCO3+   CaSO4  CaCl+  CaCl2  MgSO4   NaSO4-  H+  HCO3-  Ca+2    SO4-2    Cl-    Mg+2  Na+
    {    -1,     0,     0,      0,     0,     0,     0,     0,     0,     -1,     0,     0,     0,     0,     0,     0  }, //     OH- + H+ = H2O         
    {     0,    -1,     0,      0,     0,     0,     0,     0,     0,      1,     1,     0,     0,     0,     0,     0  }, //    CO2 + H2O = H+ + HCO3-  
    {     0,     0,    -1,      0,     0,     0,     0,     0,     0,     -1,     1,     0,     0,     0,     0,     0  }, //   CO3-2 + H+ = HCO3-       
    {     0,     0,     0,     -1,     0,     0,     0,     0,     0,      0,     1,     1,     0,     0,     0,     0  }, //      CaHCO3+ = Ca+2 + HCO3-
    {     0,     0,     0,      0,    -1,     0,     0,     0,     0,      0,     0,     1,     1,     0,     0,     0  }, //        CaSO4 = Ca+2 + SO4-2
    {     0,     0,     0,      0,     0,    -1,     0,     0,     0,      0,     0,     1,     0,     1,     0,     0  }, //        CaCl+ = Ca+2 + Cl-  
    {     0,     0,     0,      0,     0,     0,    -1,     0,     0,      0,     0,     1,     0,     2,     0,     0  }, //        CaCl2 = Ca+2 + 2Cl- 
    {     0,     0,     0,      0,     0,     0,     0,    -1,     0,      0,     0,     0,     1,     0,     1,     0  }, //        MgSO4 = Mg+2 + SO4-2
    {     0,     0,     0,      0,     0,     0,     0,     0,    -1,      0,     0,     0,     1,     0,     0,     1  }, //       NaSO4- = Na+ + SO4-2
    {     0,     0,     0,      0,     0,     0,     0,     0,     0,     -1,     1,     1,     0,     0,     0,     0  }  //   CaCO3(s) + H+ = Ca+2 + HCO3- (kinetic)
  };

// thermodynamic constants derived from 'llnl.tdat' used by Geochemists' Workbench (originally from EQ36)
constexpr CArrayWrapper<double, 10> equilibriumConstants = 
  { 
    9.89E+13,  //   OH- + H+ = H2O         
    4.42E-07,  //  CO2 + H2O = H+ + HCO3-  
    2.21E+10,  // CO3-2 + H+ = HCO3-       
    6.00E-02,  //    CaHCO3+ = Ca+2 + HCO3-
    4.79E-03,  //      CaSO4 = Ca+2 + SO4-2
    2.00E-01,  //      CaCl+ = Ca+2 + Cl-  
    3.98E+00,  //      CaCl2 = Ca+2 + 2Cl- 
    5.92E-03,  //      MgSO4 = Mg+2 + SO4-2
    2.02E-01,  //     NaSO4- = Na+ + SO4-2 
    5.16E+01   // CaCO3 + H+ = Ca+2 + HCO3- (kinetic) 
  };

constexpr CArrayWrapper<double, 10> forwardRates = 
  { 
    1.4e11,   //   OH- + H+ = H2O         
    0.039,    //  CO2 + H2O = H+ + HCO3-  
    1.0e10,   // CO3-2 + H+ = HCO3-        
    1.5e6,    //    CaHCO3+ = Ca+2 + HCO3-
    1.0e5,    //      CaSO4 = Ca+2 + SO4-2
    1.0e8,    //      CaCl+ = Ca+2 + Cl-  
    1.0e7,    //      CaCl2 = Ca+2 + 2Cl- 
    1.0e5,    //      MgSO4 = Mg+2 + SO4-2
    1.0e7,    //     NaSO4- = Na+ + SO4-2 
    1.55E-06  // CaCO3 + H+ = Ca+2 + HCO3- (kinetic) 
  };

constexpr CArrayWrapper<double, 10> reverseRates = 
  { 1.43E-03,  //   OH- + H+ = H2O         
    8.92E+04,  //  CO2 + H2O = H+ + HCO3-  
    4.67E-01,  // CO3-2 + H+ = HCO3-       
    1.85E+07,  //    CaHCO3+ = Ca+2 + HCO3-
    1.45E+07,  //      CaSO4 = Ca+2 + SO4-2
    2.14E+07,  //      CaCl+ = Ca+2 + Cl-  
    2.51E+06,  //      CaCl2 = Ca+2 + 2Cl- 
    2.69E+07,  //      MgSO4 = Mg+2 + SO4-2
    6.62E+07,  //     NaSO4- = Na+ + SO4-2
    3.00E-08   // CaCO3 + H+ = Ca+2 + HCO3-
  };

constexpr CArrayWrapper<int, 10> mobileSpeciesFlag = 
  { 1,   //   OH- + H+ = H2O         
    1,   //  CO2 + H2O = H+ + HCO3-  
    1,   // CO3-2 + H+ = HCO3-       
    1,   //    CaHCO3+ = Ca+2 + HCO3-
    1,   //      CaSO4 = Ca+2 + SO4-2
    1,   //      CaCl+ = Ca+2 + Cl-  
    1,   //      CaCl2 = Ca+2 + 2Cl- 
    1,   //      MgSO4 = Mg+2 + SO4-2
    1,   //     NaSO4- = Na+ + SO4-2
    1   // CaCO3 + H+ = Ca+2 + HCO3-
  };

// H2O coefficient, product-positive like the rows of stoichMatrix.
constexpr CArrayWrapper<signed char, 10> waterStoichiometry =
  { 1,   //   OH- + H+ = H2O
    -1,  //  CO2 + H2O = H+ + HCO3-
    0,   // CO3-2 + H+ = HCO3-
    0,   //    CaHCO3+ = Ca+2 + HCO3-
    0,   //      CaSO4 = Ca+2 + SO4-2
    0,   //      CaCl+ = Ca+2 + Cl-
    0,   //      CaCl2 = Ca+2 + 2Cl-
    0,   //      MgSO4 = Mg+2 + SO4-2
    0,   //     NaSO4- = Na+ + SO4-2
    0   // CaCO3 + H+ = Ca+2 + HCO3-
  };

// Activity model parameters
constexpr CArrayWrapper<double, 17> speciesCharge =
  // OH-      CO2(aq)   CO3-2     CaHCO3+   CaSO4(aq) CaCl+     CaCl2(aq) MgSO4(aq) NaSO4-    CaCO3(aq) H+        HCO3-     Ca+2      SO4-2     Cl-       Mg+2      Na+
  {  -1.0,    0.0,      -2.0,     1.0,      0.0,      1.0,      0.0,      0.0,      -1.0,     0.0,      1.0,      -1.0,     2.0,      -2.0,     -1.0,     2.0,      1.0 };

  // ion size parameter in ANGSTROM (phreeqc.dat -gamma values; 0.0 for neutral species and
  // species without a -gamma entry, where gamma ≈ 1)
  constexpr CArrayWrapper<double, 17> ionSize =
  // OH-      CO2(aq)   CO3-2     CaHCO3+   CaSO4(aq) CaCl+     CaCl2(aq) MgSO4(aq) NaSO4-    CaCO3(aq) H+        HCO3-     Ca+2      SO4-2     Cl-       Mg+2      Na+
  {  3.5,     0.0,      5.4,      5.4,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      9.0,      5.4,      5.0,      5.0,      3.5,      5.5,      4.0 };

  constexpr CArrayWrapper<double, 17> bdotParameters =
  // OH-      CO2(aq)   CO3-2     CaHCO3+   CaSO4(aq) CaCl+     CaCl2(aq) MgSO4(aq) NaSO4-    CaCO3(aq) H+        HCO3-     Ca+2      SO4-2     Cl-       Mg+2      Na+
  {  0.0,     0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.165,    -0.040,   0.015,    0.200,    0.075 };


// EQ3/6 B-dot parameters (data0.com.V8.R6), for validation against EQ3NR. The WATEQ form reduces
// to EQ3/6 B-dot when all species share one b. EQ3/6 applies b to charged species only.
// CO2(aq) will not match: EQ3/6 gives it a Drummond salting-out term rather than gamma = 1.
constexpr CArrayWrapper<double, 17> ionSizeEQ36 =
  // OH-      CO2(aq)   CO3-2     CaHCO3+   CaSO4(aq) CaCl+     CaCl2(aq) MgSO4(aq) NaSO4-    CaCO3(aq) H+        HCO3-     Ca+2      SO4-2     Cl-       Mg+2      Na+
  {  3.5,     3.0,      4.5,      4.0,      3.0,      4.0,      3.0,      3.0,      4.0,      3.0,      9.0,      4.0,      6.0,      4.0,      3.0,      8.0,      4.0 };

constexpr double bdotEQ36_25C = 0.0410;

constexpr CArrayWrapper<double, 17> bdotParametersEQ36 =
  // OH-           CO2(aq)       CO3-2         CaHCO3+       CaSO4(aq)     CaCl+         CaCl2(aq)     MgSO4(aq)     NaSO4-        CaCO3(aq)     H+            HCO3-         Ca+2          SO4-2         Cl-           Mg+2          Na+
  {  bdotEQ36_25C, 0.0,          bdotEQ36_25C, bdotEQ36_25C, 0.0,          bdotEQ36_25C, 0.0,          0.0,          bdotEQ36_25C, 0.0,          bdotEQ36_25C, bdotEQ36_25C, bdotEQ36_25C, bdotEQ36_25C, bdotEQ36_25C, bdotEQ36_25C, bdotEQ36_25C };

// EQ3/6 'neutral ion type' column, transcribed from the same database. 0 is neutralSpeciesType::standard,
// -1 is neutralSpeciesType::drummond.
constexpr CArrayWrapper<signed char, 17> neutralSpeciesTypeEQ36 =
  // OH-      CO2(aq)   CO3-2     CaHCO3+   CaSO4(aq) CaCl+     CaCl2(aq) MgSO4(aq) NaSO4-    CaCO3(aq) H+        HCO3-     Ca+2      SO4-2     Cl-       Mg+2      Na+
  {  0,       -1,       0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0,        0 };

}

using carbonateSystemAllKineticType     = reactionsSystems::MixedReactionsParameters< double, int, signed char, 17, 10, 0 >;
using carbonateSystemAllEquilibriumType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 17, 10, 10 >;
using carbonateSystemType               = reactionsSystems::MixedReactionsParameters< double, int, signed char, 16, 10, 9 >;

// The species count of an activity model must match that of the system it is applied to, so it is
// taken from the system type rather than repeated as a literal.
using carbonateIonicStrengthType = SpeciatedIonicStrength< double, int, carbonateSystemAllKineticType::numSpecies() >;
using carbonateActivityType = Bdot< double, int, carbonateIonicStrengthType >;
using carbonateIdentityActivityType = Identity< double, int, carbonateIonicStrengthType >;
using carbonateNosolidIonicStrengthType = SpeciatedIonicStrength< double, int, carbonateSystemType::numSpecies() >;
using carbonateNosolidActivityType = Bdot< double, int, carbonateNosolidIonicStrengthType >;
using carbonateNosolidIdentityActivityType = Identity< double, int, carbonateNosolidIonicStrengthType >;


constexpr carbonateSystemAllKineticType carbonateSystemAllKinetic( carbonate::stoichMatrix, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates, carbonate::mobileSpeciesFlag, reactionsSystems::ReactionRateLawOption::Elementary, carbonate::waterStoichiometry );
constexpr carbonateSystemAllEquilibriumType carbonateSystemAllEquilibrium( carbonate::stoichMatrix, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates, carbonate::mobileSpeciesFlag, reactionsSystems::ReactionRateLawOption::Affinity, carbonate::waterStoichiometry );
constexpr carbonateSystemType carbonateSystem( carbonate::stoichMatrixNosolid, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates, carbonate::mobileSpeciesFlag, reactionsSystems::ReactionRateLawOption::Affinity, carbonate::waterStoichiometry );

constexpr CArrayWrapper< double, 16 > carbonateNosolidSpeciesCharge =
{
  carbonate::speciesCharge[0], carbonate::speciesCharge[1], carbonate::speciesCharge[2], carbonate::speciesCharge[3],
  carbonate::speciesCharge[4], carbonate::speciesCharge[5], carbonate::speciesCharge[6], carbonate::speciesCharge[7],
  carbonate::speciesCharge[8], carbonate::speciesCharge[10], carbonate::speciesCharge[11], carbonate::speciesCharge[12],
  carbonate::speciesCharge[13], carbonate::speciesCharge[14], carbonate::speciesCharge[15], carbonate::speciesCharge[16]
};

// ion size parameter in ANGSTROM
constexpr CArrayWrapper< double, 16 > carbonateNosolidIonSize =
{
  carbonate::ionSize[0], carbonate::ionSize[1], carbonate::ionSize[2], carbonate::ionSize[3],
  carbonate::ionSize[4], carbonate::ionSize[5], carbonate::ionSize[6], carbonate::ionSize[7],
  carbonate::ionSize[8], carbonate::ionSize[10], carbonate::ionSize[11], carbonate::ionSize[12],
  carbonate::ionSize[13], carbonate::ionSize[14], carbonate::ionSize[15], carbonate::ionSize[16]
};

constexpr CArrayWrapper< double, 16 > carbonateNosolidBdotParameters =
{
  carbonate::bdotParameters[0], carbonate::bdotParameters[1], carbonate::bdotParameters[2], carbonate::bdotParameters[3],
  carbonate::bdotParameters[4], carbonate::bdotParameters[5], carbonate::bdotParameters[6], carbonate::bdotParameters[7],
  carbonate::bdotParameters[8], carbonate::bdotParameters[10], carbonate::bdotParameters[11], carbonate::bdotParameters[12],
  carbonate::bdotParameters[13], carbonate::bdotParameters[14], carbonate::bdotParameters[15], carbonate::bdotParameters[16]
};

constexpr carbonateActivityType::Params carbonateActivityParams = 
{
  {carbonate::speciesCharge},
  carbonate::ionSize,
  carbonate::bdotParameters
};

constexpr carbonateNosolidActivityType::Params carbonateNosolidActivityParams =
{
  {carbonateNosolidSpeciesCharge},
  carbonateNosolidIonSize,
  carbonateNosolidBdotParameters
};

constexpr CArrayWrapper< double, 16 > carbonateNosolidIonSizeEQ36 =
{
  carbonate::ionSizeEQ36[0], carbonate::ionSizeEQ36[1], carbonate::ionSizeEQ36[2], carbonate::ionSizeEQ36[3],
  carbonate::ionSizeEQ36[4], carbonate::ionSizeEQ36[5], carbonate::ionSizeEQ36[6], carbonate::ionSizeEQ36[7],
  carbonate::ionSizeEQ36[8], carbonate::ionSizeEQ36[10], carbonate::ionSizeEQ36[11], carbonate::ionSizeEQ36[12],
  carbonate::ionSizeEQ36[13], carbonate::ionSizeEQ36[14], carbonate::ionSizeEQ36[15], carbonate::ionSizeEQ36[16]
};

constexpr CArrayWrapper< double, 16 > carbonateNosolidBdotParametersEQ36 =
{
  carbonate::bdotParametersEQ36[0], carbonate::bdotParametersEQ36[1], carbonate::bdotParametersEQ36[2], carbonate::bdotParametersEQ36[3],
  carbonate::bdotParametersEQ36[4], carbonate::bdotParametersEQ36[5], carbonate::bdotParametersEQ36[6], carbonate::bdotParametersEQ36[7],
  carbonate::bdotParametersEQ36[8], carbonate::bdotParametersEQ36[10], carbonate::bdotParametersEQ36[11], carbonate::bdotParametersEQ36[12],
  carbonate::bdotParametersEQ36[13], carbonate::bdotParametersEQ36[14], carbonate::bdotParametersEQ36[15], carbonate::bdotParametersEQ36[16]
};

constexpr CArrayWrapper< signed char, 16 > carbonateNosolidNeutralSpeciesTypeEQ36 =
{
  carbonate::neutralSpeciesTypeEQ36[0], carbonate::neutralSpeciesTypeEQ36[1], carbonate::neutralSpeciesTypeEQ36[2], carbonate::neutralSpeciesTypeEQ36[3],
  carbonate::neutralSpeciesTypeEQ36[4], carbonate::neutralSpeciesTypeEQ36[5], carbonate::neutralSpeciesTypeEQ36[6], carbonate::neutralSpeciesTypeEQ36[7],
  carbonate::neutralSpeciesTypeEQ36[8], carbonate::neutralSpeciesTypeEQ36[10], carbonate::neutralSpeciesTypeEQ36[11], carbonate::neutralSpeciesTypeEQ36[12],
  carbonate::neutralSpeciesTypeEQ36[13], carbonate::neutralSpeciesTypeEQ36[14], carbonate::neutralSpeciesTypeEQ36[15], carbonate::neutralSpeciesTypeEQ36[16]
};

constexpr carbonateActivityType::Params carbonateActivityParamsEQ36 =
{
  {carbonate::speciesCharge},
  carbonate::ionSizeEQ36,
  carbonate::bdotParametersEQ36,
  carbonate::neutralSpeciesTypeEQ36,
  carbonate::bdotEQ36_25C // the single b the water activity assumes all solutes share
};

constexpr carbonateNosolidActivityType::Params carbonateNosolidActivityParamsEQ36 =
{
  {carbonateNosolidSpeciesCharge},
  carbonateNosolidIonSizeEQ36,
  carbonateNosolidBdotParametersEQ36,
  carbonateNosolidNeutralSpeciesTypeEQ36,
  carbonate::bdotEQ36_25C // the single b the water activity assumes all solutes share
};



constexpr Identity< double, int, carbonateIonicStrengthType >::Params carbonateIdentityActivityParams = {};
constexpr Identity< double, int, carbonateNosolidIonicStrengthType >::Params carbonateNosolidIdentityActivityParams = {};

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact
