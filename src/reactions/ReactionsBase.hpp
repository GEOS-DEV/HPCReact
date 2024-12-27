
#ifndef REACTIONS_REACTIONSBASE_HPP
#define REACTIONS_REACTIONSBASE_HPP

#include <cmath>

namespace hpcReact
{




template< typename REAL_TYPE, 
          template< int > typename REAL_DATA_ARRAY_VIEW_TYPE,
          template< int > typename REAL_CONST_DATA_ARRAY_VIEW_TYPE,
          typename INT_TYPE,
          typename INT_DATA_ARRAY_VIEW_TYPE,
          typename INT_CONST_DATA_ARRAY_VIEW_TYPE,
          typename INDEX_TYPE
        >
class ReactionsBase
{
public:
  using RealType = REAL_TYPE;
  using RealDataArrayView1d = REAL_DATA_ARRAY_VIEW_TYPE<1>;
  using RealConstDataArrayView1d = REAL_CONST_DATA_ARRAY_VIEW_TYPE<1>;
  using RealDataArrayView2d = REAL_DATA_ARRAY_VIEW_TYPE<2>;
  using RealConstDataArrayView2d = REAL_CONST_DATA_ARRAY_VIEW_TYPE<2>;
  using IntType = INT_TYPE;
  using IntDataArrayViewType = INT_DATA_ARRAY_VIEW_TYPE;
  using IntConstDataArrayType = INT_CONST_DATA_ARRAY_VIEW_TYPE;
  using IndexType = INDEX_TYPE;


  template< typename PARAMS_DATA >
  static void computeLog10ActCoefBDotModel( REAL_TYPE const temperature,
                                     REAL_TYPE const ionicStrength,
                                     PARAMS_DATA const & params,
                                     RealDataArrayView1d & log10PrimaryActCoeff,
                                     RealDataArrayView1d & dLog10PrimaryActCoeff_dIonicStrength,
                                     RealDataArrayView1d & log10SecActCoeff,
                                     RealDataArrayView1d & dLog10SecActCoeff_dIonicStrength );

  template< typename PARAMS_DATA >
  static void computeIonicStrength( PARAMS_DATA const & params,
                             RealConstDataArrayView1d const & primarySpeciesConcentration,
                             RealConstDataArrayView1d const & secondarySpeciesConcentration,
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
