#pragma once

#include "Parameters.hpp"

namespace hpcReact
{
namespace bulkGeneric
{

constexpr
KineticParameters< double, int, 5, 5 > bicarbonateBuffer =
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
  { { {  1,  1, -1,  0,  0 },
      {  1,  0, -1,  0, -1 },
      {  0, -1, -1,  1,  0 },
      {  0,  0, -1,  1,  1 },
      {  0,  1,  0,  0,  1 } }
  },
  // m_rateConstant
  { 8.42E+03, 3.71E-02, 1.25E+06, 1.24E+12, 2.30E+10 }
};

} // namespace bulkGeneric
} // namespace hpcReact
