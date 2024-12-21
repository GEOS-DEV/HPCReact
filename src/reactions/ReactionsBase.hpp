
#ifndef REACTIONS_REACTIONSBASE_HPP
#define REACTIONS_REACTIONSBASE_HPP



namespace hpcReact
{




template< typename REAL_DATA_ARRAY_VIEW_TYPE,
          typename INTEGER_DATA_ARRAY_VIEW_TYPE,
          typename REAL_TYPE, 
          typename INT_TYPE,
          typename INDEX_TYPE
        >
class ReactionsBase
{
public:
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
                             REAL_DATA_ARRAY_VIEW_TYPE const & primarySpeciesConcentration,
                             REAL_DATA_ARRAY_VIEW_TYPE const & secondarySpeciesConcentration,
                             REAL_TYPE & ionicStrength );


  template< int NUM_PRIMARY_SPECIES, 
            int NUM_SECONDARY_SPECIES >
  struct ParamsData
  {
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
