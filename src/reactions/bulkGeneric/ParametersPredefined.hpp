#pragma once

#include "Parameters.hpp"

namespace hpcReact
{
namespace bulkGeneric
{

constexpr
KineticParameters< double, int, int, 6, 5 > bicarbonateBuffer =
{
  // taken from
  // Simplified Models of the Bicarbonate Buffer for Scaled Simulations of CO2 Electrolyzers
  // Thomas Moore, Tiras Y. Lin, Thomas Roy, Sarah E. Baker, Eric B. Duoss, Christopher Hahn, and Victor A. Beck
  // Industrial & Engineering Chemistry Research 2023 62 (40), 16291-16301
  // DOI: 10.1021/acs.iecr.3c02504

  // C02, OH-, HCO3-, CO3^2-, H+
  // m_activationEnergy
  { 0.0, 0.0, 0.0, 0.0, 0.0 },
  // m_equilibriumConstant
  { 4.27E+07, 4.27E-07, 2.09E-04, 2.09E+10, 1.00E+14 },
  // m_stoichiometricMatrix
  { { -1, -1, 1, 0, 0, 0 },
    { -1, 0, 1, 0, 1, -1 },
    {  0, 1, 1, -1, 0, 0 },
    {  0, 0, 1, -1, -1, 0 },
    {  0, -1, 0, 0, -1, 1 }
  },
  // m_rateConstant
  { 8.42E+03, 3.71E-02, 1.25E+06, 1.24E+12, 2.30E+10 }
};


constexpr
KineticParameters< double, int, int, 6, 5 > bicarbonateBufferTest =
{
  // taken from
  // Simplified Models of the Bicarbonate Buffer for Scaled Simulations of CO2 Electrolyzers
  // Thomas Moore, Tiras Y. Lin, Thomas Roy, Sarah E. Baker, Eric B. Duoss, Christopher Hahn, and Victor A. Beck
  // Industrial & Engineering Chemistry Research 2023 62 (40), 16291-16301
  // DOI: 10.1021/acs.iecr.3c02504

  // m_activationEnergy
  { 0.0, 0.0, 0.0, 0.0, 0.0 },
  // m_equilibriumConstant
  { 4.27E+07, 4.27E-07, 2.09E-04, 2.09E+10, 1.00E+14 },
  // m_stoichiometricMatrix
  // C02, OH-, HCO3-, CO3^2-, H+, H20
  { { -1, -1, 1, 0, 0, 0 },
    { -1, 0, 1, 0, 1, -1 },
    {  0, 1, 1, -1, 0, 0 },
    {  0, 0, 1, -1, -1, 0 },
    {  0, -1, 0, 0, -1, 1 }
  },
  // m_rateConstant
  { 8.42E+03, 3.71E-02, 1.25E+03, 1.24E+3, 2.30E+3 }
};

constexpr
KineticParameters< double, int, int, 5, 2 > simpleTest =
{
  // m_activationEnergy
  { 0.0, 0.0 },
  // m_equilibriumConstant
  { 1.0, 1.0 },
  // m_stoichiometricMatrix
  // C02, OH-, HCO3-, CO3^2-, H+, H20
  { { -2, 1, 1, 0, 0 },
    {  0, 0, -1, -1, 2 }
  },
  // m_rateConstant
  { 1.0, 0.5 }
};

constexpr
KineticParameters< double, int, int, 10, 5 > um0 =
{
  // R1: C02 + H2O             <-> H2CO3
  // R2: H2CO3                 <-> H+ + HCO3-
  // R3: HCO3-                 <-> H+ + CO3^2-
  // R4: Mg2SiO4(solid) + 4 H+ <-> 2 Mg2+ + SiO2 + 2 H2O
  // R5: Mg^2+ + CO3^2-        <-> MgCO3(solid)

  // include Fe? no
  // add serpentine option

  //R1 = - k1f[CO2][H2O] + k1r[H2CO3] = 0   ==> H2CO3 = K1 [CO2][H2O]
  //R2 = - k2f[H2CO3] + k2r[H+][HCO3-] = 0  ==> - K2 [H2CO3] = [H+][HCO3-] ==> - K1 K2 [CO2][H2O] = [H+][HCO3-]
  //R3 = - k3f[HCO3-] + k3r[H+][CO3^2-] = 0 ==> - K3 [HCO3-] = [H+][CO3^2-] ==>

  // m_activationEnergy
//    { 79.0e3, 50.0e3, 50.0e3, 75.0e3, 30.0e3 },
  { 0, 0, 0, 0, 0 },
  // m_equilibriumConstant
  { 1.70e-3, 3.2e-4, 4.7e-11, 1.0e7, 1.0e4 },

  // m_stoichiometricMatrix
  // C02, H2CO3, H+, HCO3-, CO3^2-, Mg2SiO4, Mg2+, SiO2, MgCO3, H20
  { { -1, 1, 0, 0, 0, 0, 0, 0, 0, -1 },
    {  0, -1, 1, 1, 0, 0, 0, 0, 0, 0 },
    {  0, 0, 1, -1, 1, 0, 0, 0, 0, 0 },
    {  0, 0, -4, 0, 0, -1, 2, 1, 0, 2 },
    {  0, 0, 0, 0, -1, 0, -1, 0, 1, 0 }
  },
  // m_rateConstant
  { 0.039, 1.0e4, 1.0e4, 1.0e-12, 1.0e8 }

};


constexpr
EquilibriumKineticsModelConstants< double, int, int, 5 >
um1Constants =
{
  { 1.70e-3, 3.2e-4, 4.7e-11, 1.0e7, 1.0e4 },
  {   0.039, 1.0e4, 1.0e4, 1.0e-12, 1.0e8 },
  {     0.0, 0.0, 0.0, 0.0, 0.0 },
};

constexpr
ReactionsParameters< double, int, int, EquilibriumKineticsModelConstants, 10, 5 >
um1Params =
{
  {
    {
      { -1, 1, 0, 0, 0, 0, 0, 0, 0, -1 },
      {  0, -1, 1, 1, 0, 0, 0, 0, 0, 0 },
      {  0, 0, 1, -1, 1, 0, 0, 0, 0, 0 },
      {  0, 0, -4, 0, 0, -1, 2, 1, 0, 2 },
      {  0, 0, 0, 0, -1, 0, -1, 0, 1, 0 }
    }
  },
  um1Constants };


constexpr
EquilibriumKineticsModelConstants< double, int, int, 2 >
simpleTestRateConstants =
{
  { 1.0, 1.0 },
  { 1.0, 0.5 },
  { 1.0, 0.5 },
};

constexpr
ReactionsParameters< double, int, int, EquilibriumKineticsModelConstants, 5, 2 >
simpleTestRateParams =
{
  {
    {
      { -2, 1, 1, 0, 0 },
      {  0, 0, -1, -1, 2 }
    }
  },
  simpleTestRateConstants };



} // namespace bulkGeneric
} // namespace hpcReact
