#pragma once

#include "common/constants.hpp"
#include "common/macros.hpp"

#include <cmath>

/**
 * @file DebyeHuckel.hpp
 * @brief Debye–Hückel A^γ and B parameters for aqueous electrolytes.
 *
 * This header provides helper functions to compute the Debye–Hückel
 * parameters A^γ and B in their "native" (natural-log) form for
 * molal (mol/kg) activity-coefficient models.
 *
 * The functions are expressed in terms of fundamental physical constants
 * and water properties (density and relative permittivity). They can be
 * used directly in Debye–Hückel or extended Debye–Hückel/B-dot models.
 */

/**
 * @brief Debye-Huckel A^gamma and B parameters, and the extended Debye-Huckel log10(gamma).
 * @tparam REAL_TYPE floating point type.
 */
template< typename REAL_TYPE >
class DebyeHuckel
{
public:
  /// alias for the floating point type used in the class.
  using RealType = REAL_TYPE;

  /// π (pi).
  static constexpr RealType pi        = 3.141592653589793e+00;

  /// Vacuum permittivity ε₀ [F/m].
  static constexpr RealType e0        = 8.854187812800001e-12;

  /// Elementary charge e [C].
  static constexpr RealType eChg      = 1.602176634000000e-19;

  /// Boltzmann constant k_B [J/K].
  static constexpr RealType kB        = 1.380649000000000e-23;

  /// Avogadro constant N_A [1/mol].
  static constexpr RealType NA        = 6.022140760000000e+23;


  // -------------------------------------------------------------
  // Debye–Hückel A^γ (natural log, molal scale)
  // -------------------------------------------------------------

  /**
   * @brief Debye–Hückel A^γ parameter in natural-log form.
   *
   * Computes the coefficient A^γ(T,ρ,ε_r) used in the Debye–Hückel
   * expression for the natural logarithm of the activity coefficient:
   *
   * \f[
   *  \ln \gamma_i =
   *    - A^\gamma_{\ln}(T,P) \, z_i^2
   *      \frac{\sqrt{I}}{1 + B(T,P)\, a_i \sqrt{I}}
   * \f]
   *
   * where:
   *   - \f$ I \f$ is ionic strength in mol/kg (molal),
   *   - \f$ z_i \f$ is the ionic charge,
   *   - \f$ a_i \f$ is the ion-size parameter (length),
   *   - A^γ is independent of the log base (this function is for ln).
   *
   * The implementation follows the "native" Debye–Hückel form,
   * using fundamental physical constants without any 1/ln(10) factors.
   *
   * @param T_K          Temperature in kelvin [K].
   * @param rho_w  Density of water in g/L (≈ kg/m³ numerically).
   * @param eps_r        Relative permittivity (dielectric constant) of water.
   * @return A^γ in units consistent with molal ionic strength, for use
   *         in ln(γ) expressions.
   */
  static inline HPCREACT_HOST_DEVICE
  RealType A_gamma( RealType const T_K,
                    RealType const rho_w,
                    RealType const eps_r )
  {
    RealType const num = ::pow( eChg, 3.0 ) * ::sqrt( 2.0 * pi * NA * rho_w );
    RealType const den = ::pow( 4.0 * pi * e0 * eps_r * kB * T_K, 1.5 );
    return num / den;
  }


  // -------------------------------------------------------------
  // Debye–Hückel B (natural log, molal scale)
  // -------------------------------------------------------------

  /**
   * @brief Debye–Hückel B parameter in natural-log form.
   *
   * Computes the Debye–Hückel length-scale parameter B(T,ρ,ε_r) used
   * in the extended Debye–Hückel law:
   *
   * \f[
   *  \ln \gamma_i =
   *    - A^\gamma_{\ln}(T,P) \, z_i^2
   *      \frac{\sqrt{I}}{1 + B(T,P)\, a_i \sqrt{I}} \; ,
   * \f]
   *
   * where:
   *   - \f$ I \f$ is ionic strength in mol/kg,
   *   - \f$ a_i \f$ is an ion-size parameter (length).
   *
   * The combination \f$ B a_i \sqrt{I} \f$ is dimensionless; the
   * absolute units of B therefore depend on the length units chosen
   * for \f$ a_i \f$.
   *
   * @param T_K          Temperature in kelvin [K].
   * @param rho_w  Density of water in g/L (≈ kg/m³ numerically).
   * @param eps_r        Relative permittivity (dielectric constant) of water.
   * @return B parameter for use in ln(γ) expressions.
   */
  static inline HPCREACT_HOST_DEVICE
  RealType B_gamma( RealType const T_K,
                    RealType const rho_w,
                    RealType const eps_r )
  {
    RealType const num = 2.0 * NA * rho_w * eChg * eChg;
    RealType const den = e0 * eps_r * kB * T_K;
    return ::sqrt( num / den );
  }


  /**
   * @brief Extended Debye-Huckel log10(gamma) for a single species, with A and B evaluated from
   *        the water properties.
   * @param sqrtI Square root of the molal ionic strength I.
   * @param zi Charge of the species.
   * @param ai Ion-size parameter of the species, in ANGSTROM.
   * @param T_K Temperature in kelvin [K].
   * @param rho_w Density of water in g/L (≈ kg/m³ numerically).
   * @param eps_r Relative permittivity (dielectric constant) of water.
   * @param dlog10_gamma_dI [out] d log10(gamma) / dI. Singular at I = 0.
   * @return log10(gamma) for the species.
   *
   * A_gamma() and B_gamma() return their natural-log, SI forms, so A is converted to the log10
   * scale and B is scaled to Angstrom here, to match the units of @p ai. A caller evaluating many
   * species at one temperature should instead hoist A_gamma() and B_gamma() out of its loop and
   * call the overload taking A and B directly, as Bdot does.
   */
  static inline HPCREACT_HOST_DEVICE
  RealType log10_gamma( RealType const sqrtI,
                        RealType const zi,
                        RealType const ai,
                        RealType const T_K,
                        RealType const rho_w,
                        RealType const eps_r,
                        RealType & dlog10_gamma_dI )
  {
    RealType const A = A_gamma( T_K, rho_w, eps_r ) * hpcReact::constants::invln10;
    RealType const B = B_gamma( T_K, rho_w, eps_r ) * hpcReact::constants::metersPerAngstrom;
    return log10_gamma( sqrtI, zi, ai, A, B, dlog10_gamma_dI );
  }


  /**
   * @brief Extended Debye-Huckel log10(gamma) for a single species, with A and B supplied.
   * @param sqrtI Square root of the molal ionic strength I.
   * @param zi Charge of the species.
   * @param ai Ion-size parameter of the species, in the length units of @p B.
   * @param A Debye-Huckel A parameter, on the log10 scale.
   * @param B Debye-Huckel B parameter, scaled so that B*ai*sqrt(I) is dimensionless.
   * @param dlog10_gamma_dI [out] d log10(gamma) / dI. Singular at I = 0.
   * @return log10(gamma) for the species.
   */
  static inline HPCREACT_HOST_DEVICE
  RealType log10_gamma( RealType const sqrtI,
                        RealType const zi,
                        RealType const ai,
                        RealType const A,
                        RealType const B,
                        RealType & dlog10_gamma_dI )
  {
    RealType const denom   = 1 + B * ai * sqrtI;
    dlog10_gamma_dI = -0.5 * A * zi * zi / ( sqrtI * denom * denom );
    return -A * zi * zi * sqrtI / denom;
  }

};
