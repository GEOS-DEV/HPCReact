#pragma once

#include "common/constants.hpp"
#include "common/CArrayWrapper.hpp"


namespace hpcReact
{
namespace bulkGeneric
{




template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int NUM_PRIMARY_SPECIES,
          int NUM_KINETIC_REACTIONS >
struct KineticParameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;

  static constexpr IntType numPrimarySpecies = NUM_PRIMARY_SPECIES;
  static constexpr IntType numKineticReactions = NUM_KINETIC_REACTIONS;

  RealType m_activationEnergy[numKineticReactions];
  RealType m_equilibriumConstant[numKineticReactions];
  RealType m_stoichiometricMatrix[numKineticReactions][numPrimarySpecies];

  RealType m_rateConstant[numKineticReactions];
};

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int NUM_PRIMARY_SPECIES,
          int NUM_SECONDARY_SPECIES,
          int NUM_KINETIC_REACTIONS,
          int NUM_EQUILIBRIUM_REACTIONS >
struct BulkParameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;

  static constexpr IntType numPrimarySpecies = NUM_PRIMARY_SPECIES;
  static constexpr IntType numSecondarySpecies = NUM_SECONDARY_SPECIES;
  static constexpr IntType numKineticReactions = NUM_KINETIC_REACTIONS;
  static constexpr IntType numEquilibriumReactions = NUM_EQUILIBRIUM_REACTIONS;

};



} // namespace solidStateBattery
} // namespace hpcReact
