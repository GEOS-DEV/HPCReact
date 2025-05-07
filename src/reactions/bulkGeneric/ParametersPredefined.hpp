#pragma once

#include "Parameters.hpp"

namespace hpcReact
{
namespace bulkGeneric
{
// turn off uncrustify to allow for better readability of the parameters
// *****UNCRUSTIFY-OFF******

// constexpr 
// EquilibriumKineticsModelConstants< double, int, int, 5 > 
// um1Constants = 
// {
//   { 1.70e-3, 3.2e-4, 4.7e-11,   1.0e7, 1.0e4 },
//   {   0.039,  1.0e4,   1.0e4, 1.0e-12, 1.0e8 },
//   {     0.0,    0.0,     0.0,     0.0,   0.0 },
// };

// constexpr 
// ReactionsParameters< double, int, int, EquilibriumKineticsModelConstants, 10, 5 > 
// um1Params = 
// { 
//   { 
//     { 
//       { -1,     1,  0,     0,      0,       0,    0,    0,     0,  -1 },
//       {  0,    -1,  1,     1,      0,       0,    0,    0,     0,   0 },
//       {  0,     0,  1,    -1,      1,       0,    0,    0,     0,   0 },
//       {  0,     0, -4,     0,      0,      -1,    2,    1,     0,   2 },
//       {  0,     0,  0,     0,     -1,       0,   -1,    0,     1,   0 } 
//     }
//   }, 
//   um1Constants };


constexpr 
MixedReactionsParameters< double, int, int, 5, 2 > 
simpleTestRateParams = 
{ 
  // stoichiometric matrix
  {
    { -2, 1,  1,  0, 0 },
    {  0, 0, -1, -1, 2 } 
  },
  // equilibrium constants
  { 1.0, 1.0 },
  // forward rate constants
  { 1.0, 0.5 },
  // reverse rate constants
  { 1.0, 0.5 }
};


constexpr 
MixedReactionsParameters< double, int, int, 18, 11 > 
carbonateSystem = 
{ 
  // stoichiometric matrix
  {//   OH-    CO2  CO3-2  H2CO3 CaHCO3+ CaCO3  CaSO4  CaCl+  CaCl2  MgSO4 NaSO4-      H+  HCO3-  Ca+2  SO4-2    Cl-   Mg+2    Na+
    {    -1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0  }, //     OH- + H+ = H2O         
    {     0,    -1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0,     0,     0  }, //    CO2 + H2O = H+ + HCO3-  
    {     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     0,    -1,     1,     0,     0,     0,     0,     0  }, //   CO3-2 + H+ = HCO3-       
    {     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0,     0,     0  }, //        H2CO3 = H+ + HCO3-  
    {     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0,     0  }, //      CaHCO3+ = Ca+2 + HCO3-
    {     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,    -1,     1,     1,     0,     0,     0,     0  }, //   CaCO3 + H+ = Ca+2 + HCO3-
    {     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0  }, //        CaSO4 = Ca+2 + SO4-2
    {     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     1,     0,     1,     0,     0  }, //        CaCl+ = Ca+2 + Cl-  
    {     0,     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     1,     0,     2,     0,     0  }, //        CaCl2 = Ca+2 + 2Cl- 
    {     0,     0,     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     1,     0,     1,     0  }, //        MgSO4 = Mg+2 + SO4-2
    {     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     1,     0,     0,     1  }  //       NaSO4- = Na+ + SO4-2 
  },
  // equilibrium constants
  { 9.77E+13,   //   OH- + H+ = H2O         
    4.37E-07,   //  CO2 + H2O = H+ + HCO3-  
    2.14E+10,   // CO3-2 + H+ = HCO3-       
    1.70E-04,   //      H2CO3 = H+ + HCO3-  
    8.13E-02,   //    CaHCO3+ = Ca+2 + HCO3-
    1.17E+07,   // CaCO3 + H+ = Ca+2 + HCO3-
    6.92E-03,   //      CaSO4 = Ca+2 + SO4-2
    4.68E+00,   //      CaCl+ = Ca+2 + Cl-  
    3.98E+00,   //      CaCl2 = Ca+2 + 2Cl- 
    3.72E-03,   //      MgSO4 = Mg+2 + SO4-2
    1.51E-01 }, //     NaSO4- = Na+ + SO4-2 
  // forward rate constants
    { 1.4e11,   //   OH- + H+ = H2O         
       0.039,   //  CO2 + H2O = H+ + HCO3-  
      1.0e10,   // CO3-2 + H+ = HCO3-       
        0.57,   //      H2CO3 = H+ + HCO3-  
       1.5e6,   //    CaHCO3+ = Ca+2 + HCO3-
       1.0e5,   // CaCO3 + H+ = Ca+2 + HCO3-
       1.0e5,   //      CaSO4 = Ca+2 + SO4-2
       1.0e8,   //      CaCl+ = Ca+2 + Cl-  
       1.0e7,   //      CaCl2 = Ca+2 + 2Cl- 
       1.0e5,   //      MgSO4 = Mg+2 + SO4-2
       1.0e7    //     NaSO4- = Na+ + SO4-2 
    },
  // reverse rate constants
  { 1.43E-03,
    8.92E+04,
    4.67E-01,
    3.35E+03,
    1.85E+07,
    8.55E-03,
    1.45E+07,
    2.14E+07,
    2.51E+06,
    2.69E+07,
    6.62E+07 }
};

constexpr 
MixedReactionsParameters< double, int, int, 20, 16 > 
UltramaficSimple = 
{ 
  // stoichiometric matrix
  {//  OH-    CO2(aq)    CO3--    Mg2OH+++    Mg4(OH)++++    MgOH+    Mg2CO3++    MgCO3(aq)    MgHCO3+    Mg(H3SiO4)2    MgH2SiO4    MgH3SiO4+    H2SiO4--    H3SiO4-    H4(H2SiO4)----    H6(H2SiO4)--      H+    HCO3-    Mg++    SiO2(aq)
    {   -1,         0,       0,          0,             0,       0,          0,           0,         0,             0,          0,           0,          0,         0,                0,              0,     -1,       0,      0,          0  }, //  OH- + H+ = H2O         
    {    0,        -1,       0,          0,             0,       0,          0,           0,         0,             0,          0,           0,          0,         0,                0,              0,      1,       1,      0,          0  }, //  CO2(aq) + H2O = HCO3- + H+  
    {    0,         0,      -1,          0,             0,       0,          0,           0,         0,             0,          0,           0,          0,         0,                0,              0,     -1,       1,      0,          0  }, //  CO3-- + H+ = HCO3-       
    {    0,         0,       0,         -1,             0,       0,          0,           0,         0,             0,          0,           0,          0,         0,                0,              0,     -1,       0,      2,          0  }, //  Mg2OH+++ + H+ = 2Mg++ + H2O
    {    0,         0,       0,          0,            -1,       0,          0,           0,         0,             0,          0,           0,          0,         0,                0,              0,     -4,       0,      4,          0  }, //  Mg4(OH)++++ + 4H+ = 4Mg++ + 4H2O
    {    0,         0,       0,          0,             0,      -1,          0,           0,         0,             0,          0,           0,          0,         0,                0,              0,     -1,       0,      1,          0  }, //  MgOH+ + H+ = Mg++ + H2O
    {    0,         0,       0,          0,             0,       0,         -1,           0,         0,             0,          0,           0,          0,         0,                0,              0,     -1,       1,      2,          0  }, //  Mg2CO3++ + H+ = 2Mg++ + HCO3-
    {    0,         0,       0,          0,             0,       0,          0,          -1,         0,             0,          0,           0,          0,         0,                0,              0,     -1,       1,      1,          0  }, //  MgCO3 + H+ = Mg++ + HCO3-
    {    0,         0,       0,          0,             0,       0,          0,           0,        -1,             0,          0,           0,          0,         0,                0,              0,      0,       1,      1,          0  }, //  MgHCO3+ = Mg++ + HCO3-
    {    0,         0,       0,          0,             0,       0,          0,           0,         0,            -1,          0,           0,          0,         0,                0,              0,     -2,       0,      1,          1  }, //  Mg(H3SiO4)2 + 2H+ = Mg++ + SiO2(aq) + 4H2O
    {    0,         0,       0,          0,             0,       0,          0,           0,         0,             0,         -1,           0,          0,         0,                0,              0,     -2,       0,      1,          1  }, //  MgH2SiO4 + 2H+ = Mg++ + SiO2(aq) + 2H2O
    {    0,         0,       0,          0,             0,       0,          0,           0,         0,             0,          0,          -1,          0,         0,                0,              0,     -1,       0,      1,          1  }, //  MgH3SiO4+ + H+ = Mg++ + SiO2(aq) + 2H2O
    {    0,         0,       0,          0,             0,       0,          0,           0,         0,             0,          0,           0,         -1,         0,                0,              0,     -2,       0,      0,          1  }, //  H2SiO4-- + 2H+ = SiO2(aq) + 2H2O
    {    0,         0,       0,          0,             0,       0,          0,           0,         0,             0,          0,           0,          0,        -1,                0,              0,     -1,       0,      0,          1  }, //  H3SiO4- + H+ = SiO2(aq) + 2H2O
    {    0,         0,       0,          0,             0,       0,          0,           0,         0,             0,          0,           0,          0,         0,               -1,              0,     -4,       0,      0,          4  }, //  H4(H2SiO4)---- + 4H+ = 4SiO2(aq) + 8H2O
    {    0,         0,       0,          0,             0,       0,          0,           0,         0,             0,          0,           0,          0,         0,                0,             -1,     -2,       0,      0,          4  }  //  H6(H2SiO4)-- + 2H+ = 4SiO2 + 8H2O
  },
  // equilibrium constants
  { 9.89E+13,   //  OH- + H+ = H2O         
    4.42E-07,   //  CO2(aq) + H2O = HCO3- + H+  
    2.23E+10,   //  CO3-- + H+ = HCO3-       
    2.32E+13,   //  Mg2OH+++ + H+ = 2Mg++ + H2O
    4.47E+39,   //  Mg4(OH)++++ + 4H+ = 4Mg++ + 4H2O
    6.18E+11,   //  MgOH+ + H+ = Mg++ + H2O
    7.66E+06,   //  Mg2CO3++ + H+ = 2Mg++ + HCO3-
    2.67E+07,   //  MgCO3 + H+ = Mg++ + HCO3-
    9.77E-02,   //  MgHCO3+ = Mg++ + HCO3-
    3.45E+14,   //  Mg(H3SiO4)2 + 2H+ = Mg++ + SiO2(aq) + 4H2O
    9.49E+16,   //  MgH2SiO4 + 2H+ = Mg++ + SiO2(aq) + 2H2O
    1.96E+08,   //  MgH3SiO4+ + H+ = Mg++ + SiO2(aq) + 2H2O
    8.08E+22,   //  H2SiO4-- + 2H+ = SiO2(aq) + 2H2O
    6.44E+09,   //  H3SiO4- + H+ = SiO2(aq) + 2H2O
    5.39E+35,   //  H4(H2SiO4)---- + 4H+ = 4SiO2(aq) + 8H2O
    2.72E+13    //  H6(H2SiO4)-- + 2H+ = 4SiO2 + 8H2O
  },
  // forward rate constants (filler; carried over from EQ rxns)
  { 9.89E+13,   //  OH- + H+ = H2O         
    4.42E-07,   //  CO2(aq) + H2O = HCO3- + H+  
    2.23E+10,   //  CO3-- + H+ = HCO3-       
    2.32E+13,   //  Mg2OH+++ + H+ = 2Mg++ + H2O
    4.47E+39,   //  Mg4(OH)++++ + 4H+ = 4Mg++ + 4H2O
    6.18E+11,   //  MgOH+ + H+ = Mg++ + H2O
    7.66E+06,   //  Mg2CO3++ + H+ = 2Mg++ + HCO3-
    2.67E+07,   //  MgCO3 + H+ = Mg++ + HCO3-
    9.77E-02,   //  MgHCO3+ = Mg++ + HCO3-
    3.45E+14,   //  Mg(H3SiO4)2 + 2H+ = Mg++ + SiO2(aq) + 4H2O
    9.49E+16,   //  MgH2SiO4 + 2H+ = Mg++ + SiO2(aq) + 2H2O
    1.96E+08,   //  MgH3SiO4+ + H+ = Mg++ + SiO2(aq) + 2H2O
    8.08E+22,   //  H2SiO4-- + 2H+ = SiO2(aq) + 2H2O
    6.44E+09,   //  H3SiO4- + H+ = SiO2(aq) + 2H2O
    5.39E+35,   //  H4(H2SiO4)---- + 4H+ = 4SiO2(aq) + 8H2O
    2.72E+13    //  H6(H2SiO4)-- + 2H+ = 4SiO2 + 8H2O
  },
  // reverse rate constants (filler; carried over from EQ rxns)
  { 9.89E+13,   //  OH- + H+ = H2O         
    4.42E-07,   //  CO2(aq) + H2O = HCO3- + H+  
    2.23E+10,   //  CO3-- + H+ = HCO3-       
    2.32E+13,   //  Mg2OH+++ + H+ = 2Mg++ + H2O
    4.47E+39,   //  Mg4(OH)++++ + 4H+ = 4Mg++ + 4H2O
    6.18E+11,   //  MgOH+ + H+ = Mg++ + H2O
    7.66E+06,   //  Mg2CO3++ + H+ = 2Mg++ + HCO3-
    2.67E+07,   //  MgCO3 + H+ = Mg++ + HCO3-
    9.77E-02,   //  MgHCO3+ = Mg++ + HCO3-
    3.45E+14,   //  Mg(H3SiO4)2 + 2H+ = Mg++ + SiO2(aq) + 4H2O
    9.49E+16,   //  MgH2SiO4 + 2H+ = Mg++ + SiO2(aq) + 2H2O
    1.96E+08,   //  MgH3SiO4+ + H+ = Mg++ + SiO2(aq) + 2H2O
    8.08E+22,   //  H2SiO4-- + 2H+ = SiO2(aq) + 2H2O
    6.44E+09,   //  H3SiO4- + H+ = SiO2(aq) + 2H2O
    5.39E+35,   //  H4(H2SiO4)---- + 4H+ = 4SiO2(aq) + 8H2O
    2.72E+13    //  H6(H2SiO4)-- + 2H+ = 4SiO2 + 8H2O
  },
};

// *****UNCRUSTIFY-ON******
} // namespace bulkGeneric
} // namespace hpcReact
