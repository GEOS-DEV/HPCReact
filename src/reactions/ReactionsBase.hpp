
#ifndef REACTIONS_REACTIONSBASE_HPP
#define REACTIONS_REACTIONSBASE_HPP

#include <cmath>

namespace hpcReact
{




template< typename REAL_TYPE, 
          typename REAL_DATA_ARRAY_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_VIEW_TYPE,
          typename INT_TYPE,
          typename INT_DATA_ARRAY_VIEW_TYPE,
          typename INT_CONST_DATA_ARRAY_VIEW_TYPE,
          typename INDEX_TYPE
        >
class ReactionsBase
{
public:
  using RealType = REAL_TYPE;
  using RealDataArrayViewType = REAL_DATA_ARRAY_VIEW_TYPE;
  using RealConstDataArrayViewType = REAL_CONST_DATA_ARRAY_VIEW_TYPE;
  using IntType = INT_TYPE;
  using IntDataArrayViewType = INT_DATA_ARRAY_VIEW_TYPE;
  using IntConstDataArrayViewType = INT_CONST_DATA_ARRAY_VIEW_TYPE;
  using IndexType = INDEX_TYPE;


  template< typename PARAMS_DATA >
  static void computeLog10ActCoefBDotModel( REAL_TYPE const temperature,
                                     REAL_TYPE const ionicStrength,
                                     PARAMS_DATA const & params,
                                     REAL_DATA_ARRAY_VIEW_TYPE & log10PrimaryActCoeff,
                                     REAL_DATA_ARRAY_VIEW_TYPE & dLog10PrimaryActCoeff_dIonicStrength,
                                     REAL_DATA_ARRAY_VIEW_TYPE & log10SecActCoeff,
                                     REAL_DATA_ARRAY_VIEW_TYPE & dLog10SecActCoeff_dIonicStrength );

  template< typename PARAMS_DATA >
  static void computeIonicStrength( PARAMS_DATA const & params,
                             REAL_CONST_DATA_ARRAY_VIEW_TYPE const & primarySpeciesConcentration,
                             REAL_CONST_DATA_ARRAY_VIEW_TYPE const & secondarySpeciesConcentration,
                             REAL_TYPE & ionicStrength );


  template< int NUM_PRIMARY_SPECIES, 
            int NUM_SECONDARY_SPECIES >
  struct ParamsData
  {
    using PARENT_TYPE = ReactionsBase;

    static constexpr INT_TYPE numPrimarySpecies() { return NUM_PRIMARY_SPECIES ;}

    static constexpr INT_TYPE numSecondarySpecies() { return NUM_SECONDARY_SPECIES ;}

    REAL_TYPE m_ionSizePrimary[numPrimarySpecies()];

    REAL_TYPE m_ionSizeSec[numSecondarySpecies()];

    INT_TYPE m_chargePrimary[numPrimarySpecies()];
    INT_TYPE m_chargeSec[numSecondarySpecies()];

    REAL_TYPE m_DebyeHuckelA;
    REAL_TYPE m_DebyeHuckelB;
    REAL_TYPE m_WATEQBDot;
  };

};

} // namespace hpcReact
#endif //REACTIONS_REACTIONSBASE_HPP
