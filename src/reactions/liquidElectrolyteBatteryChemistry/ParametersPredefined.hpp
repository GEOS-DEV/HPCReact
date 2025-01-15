#pragma once

#include "Parameters.hpp"

namespace hpcReact
{
namespace geochemistry
{

constexpr
Parameters< double, int, 7, 11, 2 > chemicalReactionsParams =
{
  // m_ionSizePrimary
  { 9.00, 4.00, 6.00, 4.00, 3.00, 8.00, 4.00 },
  // m_ionSizeSec
  { 3.50, 3.00, 4.50, 3.00, 4.00, 3.00, 3.00, 4.00, 3.00, 3.00, 4.00 },
  // m_chargePrimary
  { 1, -1, 2, -2, -1, 2, 1 },
  // m_chargeSec
  { -1, 0, -2, 0, 1, 0, 0, 1, 0, 0, -1 },
  // m_DebyeHuckelA
  0.5465,
  // m_DebyeHuckelB
  0.3346,
  // m_WATEQBDot
  0.0438,
  // m_eqStoichMatrix
  // First index: 0 = OH-, 1 = CO2, 2 = CO3-2, 3 = H2CO3, 4 = CaHCO3+, 5 = CaCO3, 6 = CaSO4, 7 = CaCl+, 8 = CaCl2, 9 = MgSO4, 10 = NaSO4-
  // Second index: 0 = H+, 1 = HCO3-, 2 = Ca+2, 3 = SO4-2, 4 = Cl-, 5 = Mg+2, 6 = Na+1
  { { { -1, 0, 0, 0, 0, 0, 0 },
    {  1, 1, 0, 0, 0, 0, 0 },
    { -1, 1, 0, 0, 0, 0, 0 },
    {  1, 1, 0, 0, 0, 0, 0 },
    {  0, 1, 1, 0, 0, 0, 0 },
    { -1, 1, 1, 0, 0, 0, 0 },
    {  0, 0, 1, 1, 0, 0, 0 },
    {  0, 0, 1, 1, 0, 0, 0 },
    {  0, 0, 1, 0, 2, 0, 0 },
    {  0, 0, 0, 1, 0, 1, 0 },
    {  0, 0, 0, 1, 0, 0, 1 } }
  },
  // m_eqLog10EqConst
  { 13.99, -6.36, 10.33, -3.77, -1.09, 7.07, -2.16, 0.67, 0.60, -2.43, -0.82 },
  // m_kineticStoichMatrix
  // First index: 0 = Ca(OH)2 dissolution, 1 = CaCO3 dissolution
  // Second index: 0 = H+, 1 = HCO3-, 2 = Ca+2, 3 = SO4-2, 4 = Cl-, 5 = Mg+2, 6 = Na+1
  { { { -2, 0, 1, 0, 0, 0, 0 },
    { -1, 1, 1, 0, 0, 0, 0 } }
  },
  // m_kineticLog10EqConst
  { 20.19, 1.32 },
  // m_reactionRateConstant
  { 9.95e-1, 9.95e-3 },
  // m_specificSurfaceArea
  1.0
};

} // namespace geochemistry
} // namespace hpcReact
