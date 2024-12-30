#pragma once


namespace hpcReact
{

  template< typename REAL_TYPE,
            typename INT_TYPE,
            int NUM_PRIMARY_SPECIES,
            int NUM_SECONDARY_SPECIES >
  struct EquilibriumReactionsParameterData
  {
    using RealType = REAL_TYPE;
    using IntType = INT_TYPE;

    static constexpr INT_TYPE numPrimarySpecies = NUM_PRIMARY_SPECIES;
    static constexpr INT_TYPE numSecondarySpecies = NUM_SECONDARY_SPECIES;

    REAL_TYPE m_ionSizePrimary[numPrimarySpecies];

    REAL_TYPE m_ionSizeSec[numSecondarySpecies];

    INT_TYPE m_chargePrimary[numPrimarySpecies];
    INT_TYPE m_chargeSec[numSecondarySpecies];

    REAL_TYPE m_DebyeHuckelA;
    REAL_TYPE m_DebyeHuckelB;
    REAL_TYPE m_WATEQBDot;

    REAL_TYPE m_stoichMatrix[numSecondarySpecies][numPrimarySpecies];

  };

  template< typename REAL_TYPE,
            typename INT_TYPE,
            int NUM_PRIMARY_SPECIES,
            int NUM_SECONDARY_SPECIES,
            int NUM_KINETIC_REACTIONS >
  struct KineticReactionsParameterData : EquilibriumReactionsParameterData< REAL_TYPE, INT_TYPE, NUM_PRIMARY_SPECIES, NUM_SECONDARY_SPECIES >
  {
    using Base = EquilibriumReactionsParameterData<REAL_TYPE, INT_TYPE, NUM_PRIMARY_SPECIES, NUM_SECONDARY_SPECIES >;
    using Base::numPrimarySpecies;
    using Base::numSecondarySpecies;

    static constexpr INT_TYPE numKineticReactions = NUM_KINETIC_REACTIONS;

    using Base::m_ionSizePrimary;
    using Base::m_ionSizeSec;
    using Base::m_chargePrimary;
    using Base::m_chargeSec;
    using Base::m_DebyeHuckelA;
    using Base::m_DebyeHuckelB;
    using Base::m_WATEQBDot;
    using Base::m_stoichMatrix;

    REAL_TYPE m_reactionRateConstant[NUM_KINETIC_REACTIONS];
    REAL_TYPE m_specificSurfaceArea;
  };

}