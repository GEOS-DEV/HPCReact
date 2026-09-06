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

namespace hpcReact
{

/**
 * @brief The ideal solution activity model: every activity coefficient, and the water activity,
 *        is unity.
 * @tparam REAL_TYPE floating point type.
 * @tparam INDEX_TYPE integral type used to index the species.
 * @tparam IONIC_STRENGTH_TYPE the ionic strength model, which also supplies the base of Params.
 */
template< typename REAL_TYPE,
          typename INDEX_TYPE,
          typename IONIC_STRENGTH_TYPE >
class Identity
{
public:
  /// alias for the floating point type used in the class.
  using RealType = REAL_TYPE;

  /// alias for the integral type used to index the species.
  using IndexType = INDEX_TYPE;

  /// alias for the ionic strength model used in the class.
  using IonicStrengthType = IONIC_STRENGTH_TYPE;

  /// The ideal solution requires no parameters of its own beyond those of the ionic strength model.
  struct Params : public IONIC_STRENGTH_TYPE::Params
  {};

  /**
   * @brief Ideal solution: gamma = 1 for every species, so ln(gamma) = 0 and all derivatives vanish.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @tparam ARRAY_1D The type of the array of log activity coefficients.
   * @tparam ARRAY_2D The type of the array of log activity coefficient derivatives.
   * @tparam PARAMS The type of the activity model parameters.
   * @param params activity model parameters, unused by the ideal solution
   * @param speciesConcentrations linear concentrations c_i, unused by the ideal solution
   * @param logActivityCoefficients [out] ln(gamma_i), all zero
   * @param dLogActivityCoefficients_dConcentrations [out] d ln(gamma_i) / d c_j, all zero
   */
  template< typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename ARRAY_2D,
            typename PARAMS >
  static inline HPCREACT_HOST_DEVICE
  void
  calculateLogActivityCoefficients( PARAMS const & params,
                                    ARRAY_1D_TO_CONST const & speciesConcentrations,
                                    ARRAY_1D & logActivityCoefficients,
                                    ARRAY_2D & dLogActivityCoefficients_dConcentrations )
  {
    HPCREACT_UNUSED_VAR( params );
    HPCREACT_UNUSED_VAR( speciesConcentrations );

    constexpr IndexType numSpecies = PARAMS::numSpecies();
    for( IndexType i=0; i<numSpecies; ++i )
    {
      logActivityCoefficients[i] = 0.0;
      for( IndexType j=0; j<numSpecies; ++j )
      {
        dLogActivityCoefficients_dConcentrations[i][j] = 0.0;
      }
    }
  }

  /**
   * @brief Ideal solution: the solvent is pure, so a_w = 1 and all derivatives vanish.
   * @tparam ARRAY_1D_TO_CONST The type of the array of species concentrations.
   * @tparam ARRAY_1D The type of the array of water activity derivatives.
   * @tparam PARAMS The type of the activity model parameters.
   * @param params activity model parameters, unused by the ideal solution
   * @param speciesConcentrations linear concentrations c_i, unused by the ideal solution
   * @param dLogWaterActivity_dConcentrations [out] d ln(a_w) / d c_j, all zero
   * @return ln(a_w), which is zero
   */
  template< typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D,
            typename PARAMS >
  static inline HPCREACT_HOST_DEVICE
  REAL_TYPE
  logWaterActivity( PARAMS const & params,
                    ARRAY_1D_TO_CONST const & speciesConcentrations,
                    ARRAY_1D & dLogWaterActivity_dConcentrations )
  {
    HPCREACT_UNUSED_VAR( params );
    HPCREACT_UNUSED_VAR( speciesConcentrations );

    constexpr IndexType numSpecies = PARAMS::numSpecies();
    for( IndexType j=0; j<numSpecies; ++j )
    {
      dLogWaterActivity_dConcentrations[j] = 0.0;
    }
    return 0.0;
  }

};


} // namespace hpcReact
