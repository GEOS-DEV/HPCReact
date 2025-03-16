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


// *****UNCRUSTIFY-ON******
} // namespace bulkGeneric
} // namespace hpcReact
