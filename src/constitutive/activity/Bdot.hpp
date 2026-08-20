/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: (BSD-3-Clause)
 *
 * Copyright (c) 2025- Lawrence Livermore National Security LLC
 * All rights reserved
 *
 * See top level LICENSE files for details.
 * ------------------------------------------------------------------------------------------------------------
 */
#pragma once

#include "DebyeHuckel.hpp"
#include "common/CArrayWrapper.hpp"
#include "common/constants.hpp"

namespace hpcReact
{

template< typename REAL_TYPE,
          typename INDEX_TYPE,
          typename IONIC_STRENGTH_TYPE >
class Bdot
{
public:
  using RealType = REAL_TYPE;
  using IndexType = INDEX_TYPE;



  struct Params : public IONIC_STRENGTH_TYPE::Params
  {
    /// Ion size parameter in ANGSTROM (as tabulated by phreeqc.dat).
    CArrayWrapper< RealType, IONIC_STRENGTH_TYPE::Params::numSpecies() > m_ionSizeParameter;

    /// B-dot parameter in kg/mol, so that b*I is dimensionless.
    CArrayWrapper< RealType, IONIC_STRENGTH_TYPE::Params::numSpecies() > m_bdotParameter;
  };



  /**
   * @brief Compute ln(gamma) for every species, and its derivatives wrt linear concentration.
   * @param params activity model parameters
   * @param speciesConcentrations linear concentrations c_i
   * @param logActivityCoefficients [out] ln(gamma_i)
   * @param dLogActivityCoefficients_dConcentrations [out] d ln(gamma_i) / d c_j
   *
   * The caller composes the activity as a = c * gamma. Returning gamma rather than the activity
   * keeps gamma available to callers that need to invert it (e.g. converting a secondary species'
   * activity back to a concentration for the mole balance).
   */
  template< typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static inline HPCREACT_HOST_DEVICE
  void
  calculateLogActivityCoefficients( Params const & params,
                                    ARRAY_1D_TO_CONST const & speciesConcentrations,
                                    ARRAY_1D & logActivityCoefficients,
                                    ARRAY_2D & dLogActivityCoefficients_dConcentrations )
  {

    RealType dIonicStrength_dConcentration[ Params::numSpecies() ];
    RealType const ionicStrength = IONIC_STRENGTH_TYPE::calculate( params,
                                                                   speciesConcentrations,
                                                                   dIonicStrength_dConcentration );
    RealType const sqrtI = sqrt( ionicStrength );
    RealType const rho_w = 997.0479; // kg/m3
    RealType const eps_r = 78.54; // dimensionless
    RealType const T_K = 298.15;
    RealType const A_gamma = DebyeHuckel< RealType >::A_gamma( T_K, rho_w, eps_r );
    // A_gamma is returned in its natural-log form, while the log10_gamma equation below is
    // evaluated in log10. Convert it to the log10 scale.
    RealType const A_gamma_log10 = A_gamma * constants::invln10;

    // B_gamma*sqrt(I) is an inverse Debye length in 1/m, while m_ionSizeParameter is specified
    // in Angstrom in the parameter files (e.g. Carbonate.hpp). Scale B_gamma so that the
    // B*a*sqrt(I) group is dimensionless.
    RealType const B_gamma = DebyeHuckel< RealType >::B_gamma( T_K, rho_w, eps_r ) * constants::metersPerAngstrom;
    auto const & speciesCharge = params.m_speciesCharge;
    auto const & a = params.m_ionSizeParameter;
    auto const & b = params.m_bdotParameter;

    const IndexType numSpecies = params.numSpecies();
    for( IndexType i=0; i<numSpecies; ++i )
    {
      RealType dlog10_gamma_dI;
      RealType const DebyeHuckel_term = DebyeHuckel< RealType >::log10_gamma( sqrtI,
                                                                              speciesCharge[i],
                                                                              a[i],
                                                                              A_gamma_log10,
                                                                              B_gamma,
                                                                              dlog10_gamma_dI );
      logActivityCoefficients[i] = ( DebyeHuckel_term + b[i] * ionicStrength ) * constants::ln10;

      // d ln(gamma_i)/dc_j = ln(10) * dlog10(gamma_i)/dI * dI/dc_j.
      // dlog10_gamma_dI is singular at I = 0, where the ionic strength term is dropped.
      RealType const dLogGamma_dIonicStrength =
        ionicStrength > 0.0 ?
        constants::ln10 * ( dlog10_gamma_dI + b[i] ) :
        0.0;
      for( IndexType j=0; j<numSpecies; ++j )
      {
        dLogActivityCoefficients_dConcentrations[i][j] = dLogGamma_dIonicStrength * dIonicStrength_dConcentration[j];
      }
    }
  }

};


} // namespace hpcReact
