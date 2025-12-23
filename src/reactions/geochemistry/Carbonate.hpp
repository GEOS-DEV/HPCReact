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

}

using carbonateSystemAllKineticType     = reactionsSystems::MixedReactionsParameters< double, int, signed char, 17, 10, 0 >;
using carbonateSystemAllEquilibriumType = reactionsSystems::MixedReactionsParameters< double, int, signed char, 17, 10, 10 >;
using carbonateSystemType               = reactionsSystems::MixedReactionsParameters< double, int, signed char, 16, 10, 9 >;

constexpr carbonateSystemAllKineticType carbonateSystemAllKinetic( carbonate::stoichMatrix, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates, carbonate::mobileSpeciesFlag, 0 );
constexpr carbonateSystemAllEquilibriumType carbonateSystemAllEquilibrium( carbonate::stoichMatrix, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates, carbonate::mobileSpeciesFlag );
constexpr carbonateSystemType carbonateSystem( carbonate::stoichMatrixNosolid, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates, carbonate::mobileSpeciesFlag );

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact
