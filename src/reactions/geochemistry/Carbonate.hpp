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

constexpr CArrayWrapper<double, 11, 18> stoichMatrix = 
  { //   OH-    CO2  CO3-2  H2CO3 CaHCO3+   CaSO4  CaCl+  CaCl2  MgSO4   NaSO4- CaCO3  H+  HCO3-  Ca+2    SO4-2    Cl-    Mg+2  Na+
    {    -1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0  }, //     OH- + H+ = H2O         
    {     0,    -1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0,     0,     0  }, //    CO2 + H2O = H+ + HCO3-  
    {     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     0,    -1,     1,     0,     0,     0,     0,     0  }, //   CO3-2 + H+ = HCO3-       
    {     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0,     0,     0  }, //        H2CO3 = H+ + HCO3-  
    {     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0,     0  }, //      CaHCO3+ = Ca+2 + HCO3-
    {     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0  }, //        CaSO4 = Ca+2 + SO4-2
    {     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,     1,     0,     1,     0,     0  }, //        CaCl+ = Ca+2 + Cl-  
    {     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     1,     0,     2,     0,     0  }, //        CaCl2 = Ca+2 + 2Cl- 
    {     0,     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     1,     0,     1,     0  }, //        MgSO4 = Mg+2 + SO4-2
    {     0,     0,     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     1,     0,     0,     1  }, //       NaSO4- = Na+ + SO4-2
    {     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,    -1,    -1,     1,     1,     0,     0,     0,     0  }  //   CaCO3(s) + H+ = Ca+2 + HCO3- (kinetic)
 // {     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,    -2,     0,     1,     0,     0,     0,     0  }  //   Ca(OH)2​(s) + 2H+ = Ca2+ + 2H2​O (kinetic)
  };

constexpr CArrayWrapper<double, 11, 17> stoichMatrixNosolid = 
  { //   OH-    CO2  CO3-2  H2CO3 CaHCO3+   CaSO4  CaCl+  CaCl2  MgSO4   NaSO4-  H+  HCO3-  Ca+2    SO4-2    Cl-    Mg+2  Na+
    {    -1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     -1,     0,     0,     0,     0,     0,     0  }, //     OH- + H+ = H2O         
    {     0,    -1,     0,     0,     0,     0,     0,     0,     0,     0,      1,     1,     0,     0,     0,     0,     0  }, //    CO2 + H2O = H+ + HCO3-  
    {     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     -1,     1,     0,     0,     0,     0,     0  }, //   CO3-2 + H+ = HCO3-       
    {     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,      1,     1,     0,     0,     0,     0,     0  }, //        H2CO3 = H+ + HCO3-  
    {     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,      0,     1,     1,     0,     0,     0,     0  }, //      CaHCO3+ = Ca+2 + HCO3-
    {     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,      0,     0,     1,     1,     0,     0,     0  }, //        CaSO4 = Ca+2 + SO4-2
    {     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,      0,     0,     1,     0,     1,     0,     0  }, //        CaCl+ = Ca+2 + Cl-  
    {     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,      0,     0,     1,     0,     2,     0,     0  }, //        CaCl2 = Ca+2 + 2Cl- 
    {     0,     0,     0,     0,     0,     0,     0,     0,    -1,     0,      0,     0,     0,     1,     0,     1,     0  }, //        MgSO4 = Mg+2 + SO4-2
    {     0,     0,     0,     0,     0,     0,     0,     0,     0,    -1,      0,     0,     0,     1,     0,     0,     1  }, //       NaSO4- = Na+ + SO4-2
    {     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     -1,     1,     1,     0,     0,     0,     0  }  //   CaCO3(s) + H+ = Ca+2 + HCO3- (kinetic)
 // {     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     -2,     0,     1,     0,     0,     0,     0  }  //   Ca(OH)2​(s) + 2H+ = Ca2+ + 2H2​O (kinetic)
  };

// C^{n+1} - C^n - r( C^{n+1} ) * dt = 0
constexpr CArrayWrapper<double, 11> equilibriumConstants = 
  { 
    9.77E+13,   //   OH- + H+ = H2O         
    4.37E-07,   //  CO2 + H2O = H+ + HCO3-  
    2.14E+10,   // CO3-2 + H+ = HCO3-       
    1.70E-04,   //      H2CO3 = H+ + HCO3-  
    8.13E-02,   //    CaHCO3+ = Ca+2 + HCO3-
    6.92E-03,   //      CaSO4 = Ca+2 + SO4-2
    4.68E+00,   //      CaCl+ = Ca+2 + Cl-  
    3.98E+00,   //      CaCl2 = Ca+2 + 2Cl- 
    3.72E-03,   //      MgSO4 = Mg+2 + SO4-2
    1.51E-01,   //     NaSO4- = Na+ + SO4-2 
    1.17E+07   // CaCO3 + H+ = Ca+2 + HCO3- (kinetic)
    // 1 
  };        //   Ca(OH)2​(s) + 2H+ = Ca2+ + 2H2​O (kinetic)

constexpr CArrayWrapper<double, 11> forwardRates = 
  { 
    1.4e11,   //   OH- + H+ = H2O         
    0.039,    //  CO2 + H2O = H+ + HCO3-  
    1.0e10,   // CO3-2 + H+ = HCO3-       
    0.57,     //      H2CO3 = H+ + HCO3-  
    1.5e6,    //    CaHCO3+ = Ca+2 + HCO3-
    1.0e5,    //      CaSO4 = Ca+2 + SO4-2
    1.0e8,    //      CaCl+ = Ca+2 + Cl-  
    1.0e7,    //      CaCl2 = Ca+2 + 2Cl- 
    1.0e5,    //      MgSO4 = Mg+2 + SO4-2
    1.0e7,    //     NaSO4- = Na+ + SO4-2 
    1.0e5     // CaCO3 + H+ = Ca+2 + HCO3- (kinetic)

    // 1 
  };        //   Ca(OH)2​(s) + 2H+ = Ca2+ + 2H2​O (kinetic)

constexpr CArrayWrapper<double, 11> reverseRates = 
  { 1.43E-03,   //   OH- + H+ = H2O         
    8.92E+04,   //  CO2 + H2O = H+ + HCO3-  
    4.67E-01,   // CO3-2 + H+ = HCO3-       
    3.35E+03,   //      H2CO3 = H+ + HCO3-  
    1.85E+07,   //    CaHCO3+ = Ca+2 + HCO3-
    1.45E+07,   //      CaSO4 = Ca+2 + SO4-2
    2.14E+07,   //      CaCl+ = Ca+2 + Cl-  
    2.51E+06,   //      CaCl2 = Ca+2 + 2Cl- 
    2.69E+07,   //      MgSO4 = Mg+2 + SO4-2
    6.62E+07,   //     NaSO4- = Na+ + SO4-2
    8.55E-03    // CaCO3 + H+ = Ca+2 + HCO3-
    //  1           //   Ca(OH)2​(s) + 2H+ = Ca2+ + 2H2​O (kinetic)
  };
  }
  using carbonateSystemAllKineticType     = reactionsSystems::MixedReactionsParameters< double, int, int, 18, 11, 0 >;
  using carbonateSystemAllEquilibriumType = reactionsSystems::MixedReactionsParameters< double, int, int, 18, 11, 11 >;
  using carbonateSystemType               = reactionsSystems::MixedReactionsParameters< double, int, int, 17, 11, 10 >;

  constexpr carbonateSystemAllKineticType carbonateSystemAllKinetic( carbonate::stoichMatrix, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates );
  constexpr carbonateSystemAllEquilibriumType carbonateSystemAllEquilibrium( carbonate::stoichMatrix, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates );
  constexpr carbonateSystemType carbonateSystem( carbonate::stoichMatrixNosolid, carbonate::equilibriumConstants, carbonate::forwardRates, carbonate::reverseRates );

// *****UNCRUSTIFY-ON******
} // namespace geochemistry
} // namespace hpcReact
