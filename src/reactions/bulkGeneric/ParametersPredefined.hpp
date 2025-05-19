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


using simpleKineticTestType = MixedReactionsParameters< double, int, int, 5, 2, 0 >;

constexpr 
simpleKineticTestType
simpleKineticTestRateParams = 
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

using simpleTestType = MixedReactionsParameters< double, int, int, 5, 2, 2 >;

constexpr 
simpleTestType
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

using carbonateSystemAllKineticType = MixedReactionsParameters< double, int, int, 18, 11, 0 >;
using carbonateSystemAllEquilibriumType = MixedReactionsParameters< double, int, int, 18, 11, 11 >;
using carbonateSystemType = MixedReactionsParameters< double, int, int, 18, 11, 10 >;

constexpr CArrayWrapper<double, 11, 18> stoichMatrix = 
  { //   OH-    CO2  CO3-2  H2CO3 CaHCO3+   CaSO4  CaCl+  CaCl2  MgSO4   NaSO4-  CaCO3  H+  HCO3-  Ca+2    SO4-2    Cl-    Mg+2  Na+
    {    -1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0  }, //     OH- + H+ = H2O         
    {     0,    -1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0,     0,     0  }, //    CO2 + H2O = H+ + HCO3-  
    {     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     0,    -1,     1,     0,     0,     0,     0,     0  }, //   CO3-2 + H+ = HCO3-       
    {     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0,     0,     0  }, //        H2CO3 = H+ + HCO3-  
    {     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0,     0  }, //      CaHCO3+ = Ca+2 + HCO3-
    {     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,     0,     1,     1,     0,     0,     0  }, //        CaSO4 = Ca+2 + SO4-2
    {     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     0,     1,     0,     1,     0,     0  }, //        CaCl+ = Ca+2 + Cl-  
    {     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     1,     0,     2,     0,     0  }, //        CaCl2 = Ca+2 + 2Cl- 
    {     0,     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     0,     1,     0,     1,     0  }, //        MgSO4 = Mg+2 + SO4-2
    {     0,     0,     0,     0,     0,     0,     0,     0,     0,    -1,     0,     0,     0,     0,     1,     0,     0,     1  }, //       NaSO4- = Na+ + SO4-2
    {     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,    -1,    -1,     1,     1,     0,     0,     0,     0  }  //   CaCO3 + H+ = Ca+2 + HCO3- (kinetic)
    // {     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,    -2,     0,     1,     0,     0,     0,     0  }  //   Ca(OH)2​(s) + 2H+ = Ca2+ + 2H2​O (kinetic)
  };

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

constexpr CArrayWrapper<double, 11> reverseRates = 
  { 1.43E-03,   //   OH- + H+ = H2O         
    8.92E+04,   //  CO2 + H2O = H+ + HCO3-  
    4.67E-01,   // CO3-2 + H+ = HCO3-       
    3.35E+03,   //      H2CO3 = H+ + HCO3-  
    8.55E-03,   //    CaHCO3+ = Ca+2 + HCO3-
    1.45E+07,   //      CaSO4 = Ca+2 + SO4-2
    2.14E+07,   //      CaCl+ = Ca+2 + Cl-  
    2.51E+06,   //      CaCl2 = Ca+2 + 2Cl- 
    2.69E+07,   //      MgSO4 = Mg+2 + SO4-2
    6.62E+07,   //     NaSO4- = Na+ + SO4-2
    1.85E+07    // CaCO3 + H+ = Ca+2 + HCO3-
    //  1           //   Ca(OH)2​(s) + 2H+ = Ca2+ + 2H2​O (kinetic)
  };

  carbonateSystemAllKineticType carbonateSystemAllKinetic( stoichMatrix, equilibriumConstants, forwardRates, reverseRates );
  carbonateSystemAllEquilibriumType carbonateSystemAllEquilibrium( stoichMatrix, equilibriumConstants, forwardRates, reverseRates );
  carbonateSystemType carbonateSystem( stoichMatrix, equilibriumConstants, forwardRates, reverseRates );

// *****UNCRUSTIFY-ON******
} // namespace bulkGeneric
} // namespace hpcReact
