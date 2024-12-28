#pragma once

#include "ReactionsBase.hpp"

namespace hpcReact
{

template< typename REAL_TYPE, 
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
class KineticReactions : public ReactionsBase< REAL_TYPE,
                                                   REAL_DATA_ARRAY_1D_VIEW_TYPE,
                                                   REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                                                   INT_TYPE,
                                                   INDEX_TYPE >
{
public:

  using Base = ReactionsBase< REAL_TYPE,
                              REAL_DATA_ARRAY_1D_VIEW_TYPE,
                              REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
                              INT_TYPE,
                              INDEX_TYPE >;

  using typename Base::RealType;
  using typename Base::RealDataArrayView1d;
  using typename Base::RealConstDataArrayView1d;
  using typename Base::IntType;
  using typename Base::IndexType;


  static constexpr RealType RConst = 0;//constants::gasConstant;

  template< int NUM_PRIMARY_SPECIES, 
            int NUM_SECONDARY_SPECIES,
            int NUM_KINETIC_REACTIONS >
  struct ParamsData : Base::template ParamsData< NUM_PRIMARY_SPECIES, NUM_SECONDARY_SPECIES >
  {
    using Base = typename Base::template ParamsData< NUM_PRIMARY_SPECIES, NUM_SECONDARY_SPECIES >;
    using Base::numPrimarySpecies;
    using Base::numSecondarySpecies;

    static constexpr INT_TYPE numKineticReactions = NUM_KINETIC_REACTIONS;
    using PARENT_TYPE = KineticReactions;

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

  template< typename PARAMS_DATA >
  static void computeReactionRates( RealType const & temperature,
                                    PARAMS_DATA const & params,
                                    RealConstDataArrayView1d & primarySpeciesConcentration,
                                    RealConstDataArrayView1d & secondarySpeciesConcentration,
                                    RealDataArrayView1d & reactionRates );


};

} // namespace hpcReact