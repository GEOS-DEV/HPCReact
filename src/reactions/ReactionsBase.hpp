
#ifndef REACTIONS_REACTIONSBASE_HPP
#define REACTIONS_REACTIONSBASE_HPP

#include <cmath>

namespace hpcReact
{




template< typename REAL_TYPE, 
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE
        >
class ReactionsBase
{
public:
  using RealType = REAL_TYPE;
  using RealDataArrayView1d = REAL_DATA_ARRAY_1D_VIEW_TYPE;
  using RealConstDataArrayView1d = REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE;
  using IntType = INT_TYPE;
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

};

} // namespace hpcReact
#endif //REACTIONS_REACTIONSBASE_HPP
