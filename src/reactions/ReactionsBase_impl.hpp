#include "ReactionsBase.hpp"

namespace hpcReact
{

template< typename REAL_DATA_ARRAY_VIEW_TYPE,
          typename INTEGER_DATA_ARRAY_VIEW_TYPE,
          typename REAL_TYPE, 
          typename INT_TYPE,
          typename INDEX_TYPE
        >
template< typename PARAMS_DATA >
void ReactionsBase<REAL_DATA_ARRAY_VIEW_TYPE,
                   INTEGER_DATA_ARRAY_VIEW_TYPE,
                   REAL_TYPE, 
                   INT_TYPE,
                   INDEX_TYPE>::computeLog10ActCoefBDotModel( REAL_TYPE const temperature,
                                                              REAL_TYPE const ionicStrength,
                                                              PARAMS_DATA const & params,
                                                              REAL_DATA_ARRAY_VIEW_TYPE const & log10PrimaryActCoeff,
                                                              REAL_DATA_ARRAY_VIEW_TYPE const & dLog10PrimaryActCoeff_dIonicStrength,
                                                              REAL_DATA_ARRAY_VIEW_TYPE const & log10SecActCoeff,
                                                              REAL_DATA_ARRAY_VIEW_TYPE const & dLog10SecActCoeff_dIonicStrength ) const
{
  // Compute log10(ActivityCoefficient) for basis and dependent species along with their
  // derivatives with respect to Ionic strength using the B-Dot Model
  // which is the same as the Extended Debye-Huckel model in GEOS.
  temperature;

  for( INDEX_TYPE i = 0; i < PARAMS_DATA::numPrimarySpecies(); ++i )
  {
    log10PrimaryActCoeff[i] = params.WATEQBDot * ionicStrength - params.DebyeHuckelA * params.chargePrimary[i] * params.chargePrimary[i] * sqrt( ionicStrength ) /
                              (1.0 + params.ionSizePrimary[i] * params.DebyeHuckelB * sqrt( ionicStrength ));
    dLog10PrimaryActCoeff_dIonicStrength[i] = params.WATEQBDot - params.DebyeHuckelA * params.chargePrimary[i] * params.chargePrimary[i] *
                                              (0.5 / sqrt( ionicStrength ) / (1.0 + params.ionSizePrimary[i] * params.DebyeHuckelB * sqrt( ionicStrength )) - 0.5 * params.ionSizePrimary[i] * params.DebyeHuckelB /
                                               (1.0 + params.ionSizePrimary[i] * params.DebyeHuckelB * sqrt( ionicStrength )) /
                                               (1.0 + params.ionSizePrimary[i] * params.DebyeHuckelB * sqrt( ionicStrength )));
//    log10PrimaryActCoeff[i] = 0;
//    dLog10PrimaryActCoeff_dIonicStrength[i] = 0;
  }
  for( INDEX_TYPE i = 0; i < PARAMS_DATA::numSecondarySpecies(); ++i )
  {
    log10SecActCoeff[i] = params.WATEQBDot * ionicStrength - params.DebyeHuckelA * params.chargeSec[i] * params.chargeSec[i] * sqrt( ionicStrength ) /
                          (1.0 + params.ionSizeSec[i] * params.DebyeHuckelB * sqrt( ionicStrength ));
    dLog10SecActCoeff_dIonicStrength[i] = params.WATEQBDot - params.DebyeHuckelA * params.chargeSec[i] * params.chargeSec[i] *
                                          (0.5 / sqrt( ionicStrength ) / (1.0 + params.ionSizeSec[i] * params.DebyeHuckelB * sqrt( ionicStrength )) - 0.5 * params.ionSizeSec[i] * params.DebyeHuckelB /
                                           (1.0 + params.ionSizeSec[i] * params.DebyeHuckelB * sqrt( ionicStrength )) /
                                           (1.0 + params.ionSizeSec[i] * params.DebyeHuckelB * sqrt( ionicStrength )));
//    log10SecActCoeff[i] = 0;
//    dLog10SecActCoeff_dIonicStrength[i] = 0;
  }
}


template< typename REAL_DATA_ARRAY_VIEW_TYPE,
          typename INTEGER_DATA_ARRAY_VIEW_TYPE,
          typename REAL_TYPE, 
          typename INT_TYPE,
          typename INDEX_TYPE
        >
template< typename PARAMS_DATA >
void ReactionsBase<REAL_DATA_ARRAY_VIEW_TYPE,
                   INTEGER_DATA_ARRAY_VIEW_TYPE,
                   REAL_TYPE, 
                   INT_TYPE,
                   INDEX_TYPE>::computeIonicStrength( PARAMS_DATA const & params,
                                                      REAL_DATA_ARRAY_VIEW_TYPE const & primarySpeciesConcentration,
                                                      REAL_DATA_ARRAY_VIEW_TYPE const & secondarySpeciesConcentration,
                                                      REAL_TYPE & ionicStrength ) const
{
  //get ionic strength
  ionicStrength = 0.0;
  // Primary species
  for( INDEX_TYPE i = 0; i < PARAMS_DATA::numPrimarySpecies(); ++i )
  {
    ionicStrength += 0.5 * params.chargePrimary[i] * params.chargePrimary[i] * primarySpeciesConcentration[i];
  }
  // Secondary species
  for( INDEX_TYPE j = 0; j < PARAMS_DATA::numSecondarySpecies(); ++j )
  {
    ionicStrength += 0.5 * params.chargeSec[j] * params.chargeSec[j] * secondarySpeciesConcentration[j];
  }
}

}