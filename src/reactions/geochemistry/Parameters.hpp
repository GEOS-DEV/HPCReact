#pragma once


namespace hpcReact
{


template< typename REAL_TYPE,
          typename INT_TYPE,
          int NUM_PRIMARY_SPECIES,
          int NUM_SECONDARY_SPECIES >
struct EquilibriumParameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;

  static constexpr IntType numPrimarySpecies = NUM_PRIMARY_SPECIES;
  static constexpr IntType numSecondarySpecies = NUM_SECONDARY_SPECIES;

  RealType m_ionSizePrimary[numPrimarySpecies];

  RealType m_ionSizeSec[numSecondarySpecies];

  IntType m_chargePrimary[numPrimarySpecies];
  IntType m_chargeSec[numSecondarySpecies];

  RealType m_DebyeHuckelA;
  RealType m_DebyeHuckelB;
  RealType m_WATEQBDot;

  CArrayWrapper< RealType, numSecondarySpecies, numPrimarySpecies>  m_StoichMatrix[numSecondarySpecies][numPrimarySpecies];
  RealType m_Log10EqConst[numSecondarySpecies];
};



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

  RealType m_ionSizePrimary[numPrimarySpecies];

  RealType m_ionSizeSec[numSecondarySpecies];

  IntType m_chargePrimary[numPrimarySpecies];
  IntType m_chargeSec[numSecondarySpecies];

  RealType m_DebyeHuckelA;
  RealType m_DebyeHuckelB;
  RealType m_WATEQBDot;

  CArrayWrapper< RealType, numKineticReactions, numPrimarySpecies>  m_StoichMatrix[numKineticReactions][numPrimarySpecies];
  RealType m_log10EqConst[numKineticReactions];

  RealType m_reactionRateConstant[numKineticReactions];
  RealType m_specificSurfaceArea;
};


template< typename REAL_TYPE,
          typename INT_TYPE,
          int NUM_PRIMARY_SPECIES,
          int NUM_SECONDARY_SPECIES,
          int NUM_KINETIC_REACTIONS >
struct Parameters
{
  using RealType = REAL_TYPE;
  using IntType = INT_TYPE;

  static constexpr IntType numPrimarySpecies = NUM_PRIMARY_SPECIES;
  static constexpr IntType numSecondarySpecies = NUM_SECONDARY_SPECIES;
  static constexpr IntType numKineticReactions = NUM_KINETIC_REACTIONS;

  RealType m_ionSizePrimary[numPrimarySpecies];
  RealType m_ionSizeSec[numSecondarySpecies];

  IntType m_chargePrimary[numPrimarySpecies];
  IntType m_chargeSec[numSecondarySpecies];

  RealType m_DebyeHuckelA;
  RealType m_DebyeHuckelB;
  RealType m_WATEQBDot;

  CArrayWrapper< RealType, numSecondarySpecies, numPrimarySpecies> m_eqStoichMatrix;
  RealType m_eqLog10EqConst[numSecondarySpecies];


  CArrayWrapper< RealType, numKineticReactions, numPrimarySpecies >  m_kineticStoichMatrix;
  RealType m_kineticlog10EqConst[numKineticReactions];

  RealType m_kineticReactionRateConstant[numKineticReactions];
  RealType m_kineticSpecificSurfaceArea;

  template< int ... PRIMARY_SPECIES, 
            int ... SECONDARY_SPECIES,
            int ... SPECIES_PERMUTATIONS >
  constexpr EquilibriumParameters< REAL_TYPE, INT_TYPE, NUM_PRIMARY_SPECIES,  NUM_SECONDARY_SPECIES > 
  equilibriumReactionsParams_impl( std::integer_sequence< int, PRIMARY_SPECIES... >, 
                                   std::integer_sequence< int, SECONDARY_SPECIES... >,
                                   std::integer_sequence< int, SPECIES_PERMUTATIONS... > )
  {
    return EquilibriumParameters< REAL_TYPE, INT_TYPE, NUM_PRIMARY_SPECIES,  NUM_SECONDARY_SPECIES >{
      {m_ionSizePrimary[ PRIMARY_SPECIES ]...},
      {m_ionSizeSec[ SECONDARY_SPECIES ]...},
      {m_chargePrimary[ PRIMARY_SPECIES ]...},
      {m_chargeSec[ SECONDARY_SPECIES ]...},
      m_DebyeHuckelA,
      m_DebyeHuckelB,
      m_WATEQBDot,
      m_eqStoichMatrix,
      {m_eqLog10EqConst[ SECONDARY_SPECIES ]...}
    };
  }

  constexpr EquilibriumParameters< REAL_TYPE, INT_TYPE, NUM_PRIMARY_SPECIES,  NUM_SECONDARY_SPECIES > 
  equilibriumReactions()
  {
    return equilibriumReactionsParams_impl( std::make_integer_sequence< int, NUM_PRIMARY_SPECIES >{}, 
                                            std::make_integer_sequence< int, NUM_SECONDARY_SPECIES >{},
                                            std::make_integer_sequence< int, NUM_PRIMARY_SPECIES*NUM_SECONDARY_SPECIES >{} );
  }
};


}
