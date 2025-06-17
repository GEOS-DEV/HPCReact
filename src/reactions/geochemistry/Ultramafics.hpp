#pragma once

#include "reactions/bulkGeneric/Parameters.hpp"

namespace hpcReact
{

namespace geochemistry
{
// turn off uncrustify to allow for better readability of the parameters
// *****UNCRUSTIFY-OFF******

// ################################## Ultramafic rxn set ##################################
namespace ultramafics
{

constexpr CArrayWrapper<double, 21, 25> stoichMatrix = 
{ //      OH-   CO2(aq)      CO3--      Mg2OH+++   Mg4(OH)++++   MgOH+   Mg2CO3++   MgCO3(aq)   MgHCO3+   Mg(H3SiO4)2      MgH2SiO4     MgH3SiO4+   H2SiO4--   H3SiO4-      H4(H2SiO4)----   H6(H2SiO4)--   Mg2SiO4   MgCO3   SiO2    Mg3Si2O5(OH)4     Mg(OH)2     H+      HCO3-      Mg++   SiO2(aq)  
    {    -1,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         0,       0,         0,           0,           0,        -1,       0,        0,        0   }, //  OH- + H+ = H2O          
    {     0,        -1,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         0,       0,         0,           0,           0,         1,       1,        0,        0   }, //  CO2(aq) + H2O = HCO3- + H+   
    {     0,         0,      -1,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         0,       0,         0,           0,           0,        -1,       1,        0,        0   }, //  CO3-- + H+ = HCO3-        
    {     0,         0,       0,         -1,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         0,       0,         0,           0,           0,        -1,       0,        2,        0   }, //  Mg2OH+++ + H+ = 2Mg++ + H2O 
    {     0,         0,       0,          0,         -1,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         0,       0,         0,           0,           0,        -4,       0,        4,        0   }, //  Mg4(OH)++++ + 4H+ = 4Mg++ + 4H2O 
    {     0,         0,       0,          0,          0,      -1,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         0,       0,         0,           0,           0,        -1,       0,        1,        0   }, //  MgOH+ + H+ = Mg++ + H2O 
    {     0,         0,       0,          0,          0,       0,         -1,           0,         0,             0,          0,           0,          0,         0,            0,              0,         0,       0,         0,           0,           0,        -1,       1,        2,        0   }, //  Mg2CO3++ + H+ = 2Mg++ + HCO3- 
    {     0,         0,       0,          0,          0,       0,          0,          -1,         0,             0,          0,           0,          0,         0,            0,              0,         0,       0,         0,           0,           0,        -1,       1,        1,        0   }, //  MgCO3 + H+ = Mg++ + HCO3- 
    {     0,         0,       0,          0,          0,       0,          0,           0,        -1,             0,          0,           0,          0,         0,            0,              0,         0,       0,         0,           0,           0,         0,       1,        1,        0   }, //  MgHCO3+ = Mg++ + HCO3- 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,            -1,          0,           0,          0,         0,            0,              0,         0,       0,         0,           0,           0,        -2,       0,        1,        1   }, //  Mg(H3SiO4)2 + 2H+ = Mg++ + SiO2(aq) + 4H2O 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,         -1,           0,          0,         0,            0,              0,         0,       0,         0,           0,           0,        -2,       0,        1,        1   }, //  MgH2SiO4 + 2H+ = Mg++ + SiO2(aq) + 2H2O 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,          -1,          0,         0,            0,              0,         0,       0,         0,           0,           0,        -1,       0,        1,        1   }, //  MgH3SiO4+ + H+ = Mg++ + SiO2(aq) + 2H2O 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,         -1,         0,            0,              0,         0,       0,         0,           0,           0,        -2,       0,        0,        1   }, //  H2SiO4-- + 2H+ = SiO2(aq) + 2H2O 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,        -1,            0,              0,         0,       0,         0,           0,           0,        -1,       0,        0,        1   }, //  H3SiO4- + H+ = SiO2(aq) + 2H2O 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,           -1,              0,         0,       0,         0,           0,           0,        -4,       0,        0,        4   }, //  H4(H2SiO4)---- + 4H+ = 4SiO2(aq) + 8H2O 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,             -1,         0,       0,         0,           0,           0,        -2,       0,        0,        4   }, //  H6(H2SiO4)-- + 2H+ = 4SiO2 + 8H2O 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,        -1,       0,         0,           0,           0,        -4,       0,        2,        1   }, //  Mg2SiO4 + 4H+ = 2Mg++ + SiO2(aq) + 2H2O 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         0,      -1,         0,           0,           0,        -1,       1,        1,        0   }, //  MgCO3 + H+ = Mg++ + HCO3- 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         0,       0,        -1,           0,           0,         0,       0,        0,        1   }, //  SiO2 = SiO2(aq) 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         0,       0,         0,          -1,           0,        -6,       0,        3,        2   }, //  Mg3Si2O5(OH)4 + 6H+ = 3Mg++ + 2SiO2(aq) + 5H2O 
    {     0,         0,       0,          0,          0,       0,          0,           0,         0,             0,          0,           0,          0,         0,            0,              0,         0,       0,         0,           0,          -1,        -2,       0,        1,        0   }  //  Mg(OH)2 + 2H+ = Mg++ + 2H2O 
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
    3.45E+14,   //  Mg(H3SiO4)2 + 2H+ = Mg++ + SiO2(aq) + 4H2O
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
    1.00E+10,   //  Mg(H3SiO4)2 + 2H+ = Mg++ + SiO2(aq) + 4H2O
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
    1.00E+10,   //  Mg(H3SiO4)2 + 2H+ = Mg++ + SiO2(aq) + 4H2O
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
};
   
  using ultramaficSystemAllKineticType     = bulkGeneric::MixedReactionsParameters< double, int, int, 25, 21, 0 >;
  using ultramaficSystemAllEquilibriumType = bulkGeneric::MixedReactionsParameters< double, int, int, 25, 21, 21 >;
  using ultramaficSystemType               = bulkGeneric::MixedReactionsParameters< double, int, int, 25, 21, 16 >;

  constexpr ultramaficSystemAllKineticType     ultramaficSystemAllKinetic( ultramafics::stoichMatrix, ultramafics::equilibriumConstants, ultramafics::forwardRates, ultramafics::reverseRates );
  constexpr ultramaficSystemAllEquilibriumType ultramaficSystemAllEquilibrium( ultramafics::stoichMatrix, ultramafics::equilibriumConstants, ultramafics::forwardRates, ultramafics::reverseRates );
  constexpr ultramaficSystemType               ultramaficSystem( ultramafics::stoichMatrix, ultramafics::equilibriumConstants, ultramafics::forwardRates, ultramafics::reverseRates );

} // namespace geochemistry
} // namespace hpcReact
