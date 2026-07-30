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



  template< typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D >
  static inline HPCREACT_HOST_DEVICE
  void
  calculateActivities( Params const & params,
                       ARRAY_1D_TO_CONST const & speciesConcentrations,
                       ARRAY_1D & activities,
                       ARRAY_2D & dActivities_dConcentrations )
  {

    RealType dIonicStrength_dConcentration[ Params::numSpecies() ];
    RealType const ionicStrength = IONIC_STRENGTH_TYPE::calculate( params,
                                                                   speciesConcentrations,
                                                                   dIonicStrength_dConcentration );
    RealType const sqrtI = sqrt( ionicStrength );
    RealType const rho_w = 997.0479; // kg/m3
    RealType const eps_r = 78.54;  // dimensionless
    RealType const T_K = 298.15;
    RealType const A_gamma = DebyeHuckel< RealType >::A_gamma( T_K, rho_w, eps_r );
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
                                                                              A_gamma,
                                                                              B_gamma,
                                                                              dlog10_gamma_dI );
      RealType const log10gamma = DebyeHuckel_term + b[i] * ionicStrength;
      RealType const gamma = pow( 10.0, log10gamma );

      activities[i] = speciesConcentrations[i] * gamma;
      dActivities_dConcentrations[i][i] = gamma;
      for( IndexType j=0; j<numSpecies; ++j )
      {
        dActivities_dConcentrations[i][j] += ( ( DebyeHuckel_term * 0.5 / ionicStrength ) + b[i] ) * dIonicStrength_dConcentration[j];
      }
    }
  }

private:
  // -------------------------------------------------------------
  // Physical constants (CODATA 2018-ish), SI units
  // -------------------------------------------------------------
  static constexpr double pi   = 3.14159265358979323846;
  static constexpr double e0   = 8.8541878128e-12; // vacuum permittivity [F/m]
  static constexpr double eChg = 1.602176634e-19; // elementary charge [C]
  static constexpr double kB   = 1.380649e-23; // Boltzmann constant [J/K]
  static constexpr double NA   = 6.02214076e23; // Avogadro [1/mol]
  static constexpr double ln10 = 2.3025850929940459;

};


} // namespace hpcReact
