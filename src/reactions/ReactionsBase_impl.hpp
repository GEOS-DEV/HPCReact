#include "ReactionsBase.hpp"

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
template< typename PARAMS_DATA >
void ReactionsBase< REAL_TYPE,
                    REAL_DATA_ARRAY_VIEW_TYPE,
                    REAL_CONST_DATA_ARRAY_VIEW_TYPE,
                    INT_TYPE,
                    INT_DATA_ARRAY_VIEW_TYPE,
                    INT_CONST_DATA_ARRAY_VIEW_TYPE,
                    INDEX_TYPE >::computeLog10ActCoefBDotModel( REAL_TYPE const ,//temperature,
                                                              REAL_TYPE const ionicStrength,
                                                              PARAMS_DATA const & params,
                                                              REAL_DATA_ARRAY_VIEW_TYPE & log10PrimaryActCoeff,
                                                              REAL_DATA_ARRAY_VIEW_TYPE & dLog10PrimaryActCoeff_dIonicStrength,
                                                              REAL_DATA_ARRAY_VIEW_TYPE & log10SecActCoeff,
                                                              REAL_DATA_ARRAY_VIEW_TYPE & dLog10SecActCoeff_dIonicStrength )
{
  // Compute log10(ActivityCoefficient) for basis and dependent species along with their
  // derivatives with respect to Ionic strength using the B-Dot Model
  // which is the same as the Extended Debye-Huckel model in GEOS.

  for( INDEX_TYPE i = 0; i < PARAMS_DATA::numPrimarySpecies(); ++i )
  {
    log10PrimaryActCoeff[i] = params.m_WATEQBDot * ionicStrength - params.m_DebyeHuckelA * params.m_chargePrimary[i] * params.m_chargePrimary[i] * sqrt( ionicStrength ) /
                              (1.0 + params.m_ionSizePrimary[i] * params.m_DebyeHuckelB * sqrt( ionicStrength ));
    dLog10PrimaryActCoeff_dIonicStrength[i] = params.m_WATEQBDot - params.m_DebyeHuckelA * params.m_chargePrimary[i] * params.m_chargePrimary[i] *
                                              (0.5 / sqrt( ionicStrength ) / (1.0 + params.m_ionSizePrimary[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )) - 0.5 * params.m_ionSizePrimary[i] * params.m_DebyeHuckelB /
                                               (1.0 + params.m_ionSizePrimary[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )) /
                                               (1.0 + params.m_ionSizePrimary[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )));
//    log10PrimaryActCoeff[i] = 0;
//    dLog10PrimaryActCoeff_dIonicStrength[i] = 0;
  }
  for( INDEX_TYPE i = 0; i < PARAMS_DATA::numSecondarySpecies(); ++i )
  {
    log10SecActCoeff[i] = params.m_WATEQBDot * ionicStrength - params.m_DebyeHuckelA * params.m_chargeSec[i] * params.m_chargeSec[i] * sqrt( ionicStrength ) /
                          (1.0 + params.m_ionSizeSec[i] * params.m_DebyeHuckelB * sqrt( ionicStrength ));
    dLog10SecActCoeff_dIonicStrength[i] = params.m_WATEQBDot - params.m_DebyeHuckelA * params.m_chargeSec[i] * params.m_chargeSec[i] *
                                          (0.5 / sqrt( ionicStrength ) / (1.0 + params.m_ionSizeSec[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )) - 0.5 * params.m_ionSizeSec[i] * params.m_DebyeHuckelB /
                                           (1.0 + params.m_ionSizeSec[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )) /
                                           (1.0 + params.m_ionSizeSec[i] * params.m_DebyeHuckelB * sqrt( ionicStrength )));
//    log10SecActCoeff[i] = 0;
//    dLog10SecActCoeff_dIonicStrength[i] = 0;
  }
}


template< typename REAL_TYPE, 
          typename REAL_DATA_ARRAY_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_VIEW_TYPE,
          typename INT_TYPE,
          typename INT_DATA_ARRAY_VIEW_TYPE,
          typename INT_CONST_DATA_ARRAY_VIEW_TYPE,
          typename INDEX_TYPE
        >
template< typename PARAMS_DATA >
void ReactionsBase< REAL_TYPE,
                    REAL_DATA_ARRAY_VIEW_TYPE,
                    REAL_CONST_DATA_ARRAY_VIEW_TYPE,
                    INT_TYPE,
                    INT_DATA_ARRAY_VIEW_TYPE,
                    INT_CONST_DATA_ARRAY_VIEW_TYPE,
                    INDEX_TYPE >::computeIonicStrength( PARAMS_DATA const & params,
                                                      REAL_CONST_DATA_ARRAY_VIEW_TYPE const & primarySpeciesConcentration,
                                                      REAL_CONST_DATA_ARRAY_VIEW_TYPE const & secondarySpeciesConcentration,
                                                      REAL_TYPE & ionicStrength )
{
  //get ionic strength
  ionicStrength = 0.0;
  // Primary species
  for( INDEX_TYPE i = 0; i < PARAMS_DATA::numPrimarySpecies(); ++i )
  {
    ionicStrength += 0.5 * params.m_chargePrimary[i] * params.m_chargePrimary[i] * primarySpeciesConcentration[i];
  }
  // Secondary species
  for( INDEX_TYPE j = 0; j < PARAMS_DATA::numSecondarySpecies(); ++j )
  {
    ionicStrength += 0.5 * params.m_chargeSec[j] * params.m_chargeSec[j] * secondarySpeciesConcentration[j];
  }
}

}