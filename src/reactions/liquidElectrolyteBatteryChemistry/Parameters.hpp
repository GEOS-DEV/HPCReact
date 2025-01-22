#pragma once

#include "common/constants.hpp"


namespace hpcReact
{
namespace solidStateBattery
{




template< typename REAL_TYPE,
          typename INT_TYPE,
          int NUM_PRIMARY_SPECIES,
          int NUM_SECONDARY_SPECIES,
          int NUM_KINETIC_REACTIONS >
struct KineticParameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;

  static constexpr IntType numPrimarySpecies = NUM_PRIMARY_SPECIES;
  static constexpr IntType numSecondarySpecies = NUM_SECONDARY_SPECIES;
  static constexpr IntType numKineticReactions = NUM_KINETIC_REACTIONS;

  RealType m_activationEnergy[numKineticReactions];
  RealType m_equilibriumConstant[numKineticReactions];
  CArrayWrapper< RealType, numKineticReactions, numPrimarySpecies>  m_stoichiometricMatrix[numKineticReactions][numPrimarySpecies];

  RealType m_rateConstant[numKineticReactions];
};




} // namespace solidStateBattery
} // namespace hpcReact
