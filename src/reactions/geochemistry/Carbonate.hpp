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





constexpr CArrayWrapper<double, 17> speciesCharge = 
  { -1.0, // OH-
     0.0, // CO2(aq)
    -2.0, // CO3-2
     1.0, // CaHCO3+
     0.0, // CaSO4(aq)
     1.0, // CaCl+
     0.0, // CaCl2(aq)
     0.0, // MgSO4(aq)
    -1.0, // NaSO4-
     0.0, // CaCO3(aq)
     1.0, // H+
    -1.0, // HCO3-
     2.0, // Ca+2
    -2.0, // SO4-2
    -1.0, // Cl-
     2.0, // Mg+2
     1.0  // Na+
    };

  // ion size parameter in ANGSTROM
  constexpr CArrayWrapper<double, 17> ionSize =
  {
    3.5,  // OH-      (from H2O = OH- + H+, -gamma 3.5 0.0)
    0.0,  // CO2(aq)  (neutral, no -gamma; typically gamma ≈ 1)
    5.4,  // CO3-2
    5.4,  // CaHCO3+
    0.0,  // CaSO4(aq) (neutral)
    0.0,  // CaCl+    (no -gamma in phreeqc.dat)
    0.0,  // CaCl2(aq) (neutral)
    0.0,  // MgSO4(aq) (neutral)
    0.0,  // NaSO4-   (no -gamma in phreeqc.dat)
    0.0,  // CaCO3(aq) (neutral)
    9.0,  // H+
    5.4,  // HCO3-    (from CO3-2 + H+ = HCO3-, -gamma 5.4 0.0)
    5.0,  // Ca+2
    5.0,  // SO4-2
    3.5,  // Cl-
    5.5,  // Mg+2
    4.0   // Na+ 
  };

  constexpr CArrayWrapper<double, 17> bdotParameters =
  {
      0.0,    // OH-
      0.0,    // CO2(aq)
      0.0,    // CO3-2
      0.0,    // CaHCO3+
      0.0,    // CaSO4(aq)
      0.0,    // CaCl+
      0.0,    // CaCl2(aq)
      0.0,    // MgSO4(aq)
      0.0,    // NaSO4-
      0.0,    // CaCO3(aq)
      0.0,    // H+
      0.0,    // HCO3-
    0.165,  // Ca+2
   -0.040,  // SO4-2
    0.015,  // Cl-
    0.200,  // Mg+2
    0.075   // Na+
  }; 


// EQ3/6 B-dot parameters (data0.com.V8.R6), for validation against EQ3NR. The WATEQ form reduces
// to EQ3/6 B-dot when all species share one b. EQ3/6 applies b to charged species only.
// CO2(aq) will not match: EQ3/6 gives it a Drummond salting-out term rather than gamma = 1.
constexpr CArrayWrapper<double, 17> ionSizeEQ36 =
  {
    3.5,  // OH-
    3.0,  // CO2(aq)
    4.5,  // CO3-2
    4.0,  // CaHCO3+
    3.0,  // CaSO4(aq)
    4.0,  // CaCl+
    3.0,  // CaCl2(aq)
    3.0,  // MgSO4(aq)
    4.0,  // NaSO4-
    3.0,  // CaCO3(aq)
    9.0,  // H+
    4.0,  // HCO3-
    6.0,  // Ca+2
    4.0,  // SO4-2
    3.0,  // Cl-
    8.0,  // Mg+2
    4.0   // Na+
  };

constexpr double bdotEQ36_25C = 0.0410;

constexpr CArrayWrapper<double, 17> bdotParametersEQ36 =
  {
    bdotEQ36_25C,  // OH-
    0.0,           // CO2(aq)
    bdotEQ36_25C,  // CO3-2
    bdotEQ36_25C,  // CaHCO3+
    0.0,           // CaSO4(aq)
    bdotEQ36_25C,  // CaCl+
    0.0,           // CaCl2(aq)
    0.0,           // MgSO4(aq)
    bdotEQ36_25C,  // NaSO4-
    0.0,           // CaCO3(aq)
    bdotEQ36_25C,  // H+
    bdotEQ36_25C,  // HCO3-
    bdotEQ36_25C,  // Ca+2
    bdotEQ36_25C,  // SO4-2
    bdotEQ36_25C,  // Cl-
    bdotEQ36_25C,  // Mg+2
    bdotEQ36_25C   // Na+
  };

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


constexpr carbonateSystemAllKineticType carbonateSystemAllKinetic( carbonate::stoichMatrix, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates, carbonate::mobileSpeciesFlag, 0 );
constexpr carbonateSystemAllEquilibriumType carbonateSystemAllEquilibrium( carbonate::stoichMatrix, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates, carbonate::mobileSpeciesFlag );
constexpr carbonateSystemType carbonateSystem( carbonate::stoichMatrixNosolid, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates, carbonate::mobileSpeciesFlag );

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

constexpr carbonateActivityType::Params carbonateActivityParamsEQ36 =
{
  {carbonate::speciesCharge},
  carbonate::ionSizeEQ36,
  carbonate::bdotParametersEQ36
};

constexpr carbonateNosolidActivityType::Params carbonateNosolidActivityParamsEQ36 =
{
  {carbonateNosolidSpeciesCharge},
  carbonateNosolidIonSizeEQ36,
  carbonateNosolidBdotParametersEQ36
};



Identity< double, int, carbonateIonicStrengthType >::Params carbonateIdentityActivityParams = {};
Identity< double, int, carbonateNosolidIonicStrengthType >::Params carbonateNosolidIdentityActivityParams = {};

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact
