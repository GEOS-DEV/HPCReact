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
#include "Drummond.hpp"
#include "common/CArrayWrapper.hpp"
#include "common/constants.hpp"

namespace hpcReact
{

/**
 * @brief Selects which model supplies a species' activity coefficient.
 *
 * The values are EQ3/6's own "neutral ion type" codes, as tabulated in the 'bdot parameters' block
 * of a data0 file, so a parameter file can transcribe that column without translating it. In
 * data0.com.V8.R6 exactly three of the 1769 aqueous species are tagged for salting-out --
 * CO2(aq), H2(aq) and O2(aq). Every other species carries the default, H2S(aq), N2(aq), NH3(aq)
 * and SO2(aq) among them.
 */
namespace neutralSpeciesType
{
/// The standard B-dot expression. It degenerates to gamma = 1 for a neutral species, whose
/// Debye-Huckel term vanishes with its charge.
constexpr signed char standard = 0;

/// Drummond (1981) salting-out polynomial, in place of the B-dot expression.
constexpr signed char drummond = -1;
}

/**
 * @brief The B-dot (Helgeson) activity model, with Drummond salting-out for the species tagged
 *        for it, and the B-dot-consistent water activity.
 * @tparam REAL_TYPE floating point type.
 * @tparam INDEX_TYPE integral type used to index the species.
 * @tparam IONIC_STRENGTH_TYPE the ionic strength model, which also supplies the base of Params.
 */
template< typename REAL_TYPE,
          typename INDEX_TYPE,
          typename IONIC_STRENGTH_TYPE >
class Bdot
{
public:
  /// alias for the floating point type used in the class.
  using RealType = REAL_TYPE;

  /// alias for the integral type used to index the species.
  using IndexType = INDEX_TYPE;

  /// alias for the ionic strength model used in the class.
  using IonicStrengthType = IONIC_STRENGTH_TYPE;


  /// The B-dot parameters, extending those the ionic strength model requires.
  struct Params : public IONIC_STRENGTH_TYPE::Params
  {
    /// Ion size parameter in ANGSTROM (as tabulated by phreeqc.dat).
    CArrayWrapper< RealType, IONIC_STRENGTH_TYPE::Params::numSpecies() > m_ionSizeParameter;

    /// B-dot parameter in kg/mol, so that b*I is dimensionless.
    CArrayWrapper< RealType, IONIC_STRENGTH_TYPE::Params::numSpecies() > m_bdotParameter;

    /// Per-species neutralSpeciesType tag. Defaults to all-standard, which is the behavior of a
    /// parameter file written before this member existed.
    CArrayWrapper< signed char, IONIC_STRENGTH_TYPE::Params::numSpecies() > m_neutralSpeciesType {};

    /// The single B-dot parameter the water activity assumes all solutes share. Defaults to 0.
    RealType m_bdotWater {};
  };

  /// Ambient water density [kg/m3], shared by the activity coefficients and the water activity.
  static constexpr RealType rho_w = 997.0479;

  /// Ambient relative permittivity of water [dimensionless].
  static constexpr RealType eps_r = 78.54;

  /// Ambient temperature [K].
  static constexpr RealType T_K = 298.15;



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
    auto const & neutralType = params.m_neutralSpeciesType;

    const IndexType numSpecies = params.numSpecies();
    for( IndexType i=0; i<numSpecies; ++i )
    {
      RealType dLogGamma_dIonicStrength;

      if( neutralType[i] == neutralSpeciesType::drummond )
      {
        logActivityCoefficients[i] = Drummond< RealType >::ln_gamma( ionicStrength,
                                                                     T_K,
                                                                     dLogGamma_dIonicStrength );
      }
      else
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
        dLogGamma_dIonicStrength =
          ionicStrength > 0.0 ?
          constants::ln10 * ( dlog10_gamma_dI + b[i] ) :
          0.0;
      }

      for( IndexType j=0; j<numSpecies; ++j )
      {
        dLogActivityCoefficients_dConcentrations[i][j] = dLogGamma_dIonicStrength * dIonicStrength_dConcentration[j];
      }
    }
  }

  /**
   * @brief Compute ln(a_w), the activity of the solvent, and its derivatives.
   * @param params activity model parameters
   * @param speciesConcentrations linear concentrations c_i, in molality
   * @param dLogWaterActivity_dConcentrations [out] d ln(a_w)/dc_j for every species j
   * @return ln(a_w)
   *
   * The B-dot-consistent form
   * \f[
   *   \log_{10} a_w = \frac{1}{\Omega} \left[ -\frac{\sum_i m_i}{\ln 10}
   *                 + \frac{2}{3} A^\gamma_{10} I^{3/2} \sigma( \mathring{a} B^\gamma \sqrt{I} )
   *                 - \dot{B} I^2 \right],
   *   \quad \sigma(x) = \frac{3}{x^3}\left( 1 + x - \frac{1}{1+x} - 2\ln(1+x) \right)
   * \f]
   * It is consistent with the B-dot gamma above when every solute is an ion sharing one hard core
   * diameter, one B-dot parameter and one z^2.
   */
  template< typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D >
  static inline HPCREACT_HOST_DEVICE
  RealType
  logWaterActivity( Params const & params,
                    ARRAY_1D_TO_CONST const & speciesConcentrations,
                    ARRAY_1D & dLogWaterActivity_dConcentrations )
  {
    RealType dIonicStrength_dConcentration[ Params::numSpecies() ];
    RealType const ionicStrength = IONIC_STRENGTH_TYPE::calculate( params,
                                                                   speciesConcentrations,
                                                                   dIonicStrength_dConcentration );

    RealType dLnWaterActivity_dSoluteMolality;
    RealType dLnWaterActivity_dIonicStrength;
    RealType const result = logWaterActivity_impl( params,
                                                   speciesConcentrations,
                                                   ionicStrength,
                                                   dLnWaterActivity_dSoluteMolality,
                                                   dLnWaterActivity_dIonicStrength );

    IndexType const numSpecies = params.numSpecies();
    for( IndexType j=0; j<numSpecies; ++j )
    {
      dLogWaterActivity_dConcentrations[j] = dLnWaterActivity_dSoluteMolality
                                             + dLnWaterActivity_dIonicStrength * dIonicStrength_dConcentration[j];
    }
    return result;
  }

private:

  /**
   * @brief The closed form above, returning its two partial derivatives rather than a gradient.
   * @param ionicStrength molal ionic strength I
   * @param dLnWaterActivity_dSoluteMolality [out] d ln(a_w) / d(sum_i m_i)
   * @param dLnWaterActivity_dIonicStrength [out] d ln(a_w) / dI
   * @return ln(a_w)
   */
  template< typename ARRAY_1D_TO_CONST >
  static inline HPCREACT_HOST_DEVICE
  RealType
  logWaterActivity_impl( Params const & params,
                         ARRAY_1D_TO_CONST const & speciesConcentrations,
                         RealType const ionicStrength,
                         RealType & dLnWaterActivity_dSoluteMolality,
                         RealType & dLnWaterActivity_dIonicStrength )
  {
    /// Hard core diameter in ANGSTROM, fixed for every solute.
    constexpr RealType hardCoreDiameter = 4.0;

    RealType soluteMolality = 0.0;
    IndexType const numSpecies = params.numSpecies();
    for( IndexType i=0; i<numSpecies; ++i )
    {
      soluteMolality += speciesConcentrations[i];
    }

    RealType const A_gamma_log10 = DebyeHuckel< RealType >::A_gamma( T_K, rho_w, eps_r ) * constants::invln10;
    RealType const B_gamma = DebyeHuckel< RealType >::B_gamma( T_K, rho_w, eps_r ) * constants::metersPerAngstrom;

    // I^(3/2)*sigma(k*sqrt(I)) reduces to (3/k^3)*h(x), which cancels both the I^(3/2) and the
    // 1/x^3 and so is finite at I = 0.
    RealType const k = hardCoreDiameter * B_gamma;
    RealType const x = k * sqrt( ionicStrength );
    RealType const onePlusX = 1.0 + x;
    RealType const h      = 1.0 + x - 1.0 / onePlusX - 2.0 * log( onePlusX );
    RealType const dh_dx  = 1.0 + 1.0 / ( onePlusX * onePlusX ) - 2.0 / onePlusX;

    RealType const debyeHuckelTerm = 2.0 * A_gamma_log10 * h / ( k * k * k );
    RealType const bdotTerm        = -params.m_bdotWater * ionicStrength * ionicStrength;

    // dh_dx/(k*x) is the I-derivative of the Debye-Huckel term; it tends to 0 with x.
    RealType const dTerms_dIonicStrength =
      ionicStrength > 0.0 ?
      A_gamma_log10 * dh_dx / ( k * x ) - 2.0 * params.m_bdotWater * ionicStrength :
      0.0;

    dLnWaterActivity_dSoluteMolality = -1.0 / constants::waterMolality;
    dLnWaterActivity_dIonicStrength  = constants::ln10 * dTerms_dIonicStrength / constants::waterMolality;

    return constants::ln10 * ( -soluteMolality * constants::invln10 + debyeHuckelTerm + bdotTerm )
           / constants::waterMolality;
  }

};


} // namespace hpcReact
