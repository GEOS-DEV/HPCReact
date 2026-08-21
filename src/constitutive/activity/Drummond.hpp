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

#include "common/macros.hpp"

/**
 * @file Drummond.hpp
 * @brief Drummond (1981) salting-out model for neutral aqueous species.
 */

namespace hpcReact
{

/**
 * @brief Drummond's (1981) salting-out polynomial for neutral aqueous species.
 *
 * \f[
 *   \ln \gamma_n = \left( c_1 + c_2 T + \frac{c_3}{T} \right) I
 *                - \left( c_4 + c_5 T \right) \frac{I}{I + 1}
 * \f]
 *
 * where I is the molal ionic strength and T is in kelvin. A neutral species has no Debye-Huckel
 * term, so without a model of this kind its gamma is 1; that misses the salting-out of a dissolved
 * gas, which reaches ~1.4 for CO2(aq) at I = 1.6.
 *
 * The coefficients are those of the 'cco2' block of data0.com.V8.R6, read in tabulated order. They
 * carry no species index: EQ3/6 applies this one set to every species it tags for salting-out.
 *
 * @tparam REAL_TYPE floating point type.
 */
template< typename REAL_TYPE >
class Drummond
{
public:
  using RealType = REAL_TYPE;

  /// Constant term of the linear-in-I group [dimensionless].
  static constexpr RealType c1 = -1.0312;
  /// Temperature coefficient of the linear-in-I group [1/K].
  static constexpr RealType c2 = 0.0012806;
  /// Inverse-temperature coefficient of the linear-in-I group [K].
  static constexpr RealType c3 = 255.9;
  /// Constant term of the saturating group [dimensionless].
  static constexpr RealType c4 = 0.4445;
  /// Temperature coefficient of the saturating group [1/K].
  static constexpr RealType c5 = -0.001606;

  /**
   * @brief Compute ln(gamma) for a neutral species, and its derivative wrt ionic strength.
   * @param ionicStrength molal ionic strength I
   * @param T_K temperature in kelvin
   * @param dLnGamma_dIonicStrength [out] d ln(gamma) / dI
   * @return ln(gamma)
   *
   * Unlike the Debye-Huckel term, this expression and its derivative are both finite at I = 0,
   * so no special case is needed there.
   */
  static inline HPCREACT_HOST_DEVICE
  RealType ln_gamma( RealType const ionicStrength,
                     RealType const T_K,
                     RealType & dLnGamma_dIonicStrength )
  {
    RealType const linearGroup     = c1 + c2 * T_K + c3 / T_K;
    RealType const saturatingGroup = c4 + c5 * T_K;
    RealType const onePlusI        = 1.0 + ionicStrength;

    dLnGamma_dIonicStrength = linearGroup - saturatingGroup / ( onePlusI * onePlusI );

    return linearGroup * ionicStrength - saturatingGroup * ionicStrength / onePlusI;
  }
};

} // namespace hpcReact
