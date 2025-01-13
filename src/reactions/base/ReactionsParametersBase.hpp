#pragma once


namespace hpcReact
{

template< typename REAL_TYPE,
          typename INT_TYPE,
          int NUM_PRIMARY_SPECIES,
          int NUM_SECONDARY_SPECIES >
struct ParameterBase
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;

  static constexpr IntType numPrimarySpecies = NUM_PRIMARY_SPECIES;
  static constexpr IntType numSecondarySpecies = NUM_SECONDARY_SPECIES;

  RealType m_ionSizePrimary[numPrimarySpecies];
  RealType m_ionSizeSec[numSecondarySpecies];

  IntType m_chargePrimary[numPrimarySpecies];
  IntType m_chargeSec[numSecondarySpecies];
};

}
