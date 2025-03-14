#pragma once

#include "common/macros.hpp"
#include "common/CArrayWrapper.hpp"
#include "common/DirectSystemSolve.hpp"
#include "common/printers.hpp"

#include <iostream>


namespace hpcReact
{
namespace bulkGeneric
{

/**
 * @brief This class implements components required to calculate equilibrium 
 *        reactions for a given set of species. The class also proovides
 *        device callable methods to enforce equilibrium pointwise and also
 *        provides methods to launch batch processing of material points.
 * @tparam REAL_TYPE The type of the real numbers used in the class.
 * @tparam INT_TYPE The type of the integers used in the class.
 * @tparam INDEX_TYPE The type of the indices used in the class.
 */
template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
class EquilibriumReactions
{
public:
  /// alias for type of the real numbers used in the class.
  using RealType = REAL_TYPE;

  /// alias for type of the integers used in the class.
  using IntType = INT_TYPE;
  
  /// alias for type of the indices used in the class.
  using IndexType = INDEX_TYPE;

  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST2,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeResidualAndJacobian( RealType const & temperature,
                              PARAMS_DATA const & params,
                              ARRAY_1D_TO_CONST const & speciesConcentration0,
                              ARRAY_1D_TO_CONST2 const & xi,
                              ARRAY_1D & residual,
                              ARRAY_2D & jacobian );

  /**
   * @brief This method enforces equilibrium for a given set of species using 
   *        reaction extents.
  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST >
  static HPCREACT_HOST_DEVICE
  void
  enforceEquilibrium_Extents( RealType const & temperature,
                              PARAMS_DATA const & params,
                              ARRAY_1D_TO_CONST const & speciesConcentration0,
                              ARRAY_1D & speciesConcentration );

  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST >
  static HPCREACT_HOST_DEVICE
  void
  enforceEquilibrium_Aggregate( RealType const & temperature,
                                PARAMS_DATA const & params,
                                ARRAY_1D_TO_CONST const & speciesConcentration0,
                                ARRAY_1D & speciesConcentration );

private:
  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_1D_TO_CONST2,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeResidualAndJacobianReactionExtents( RealType const & temperature,
                                             PARAMS_DATA const & params,
                                             ARRAY_1D_TO_CONST const & speciesConcentration0,
                                             ARRAY_1D_TO_CONST2 const & xi,
                                             ARRAY_1D & residual,
                                             ARRAY_2D & jacobian );

  template< typename PARAMS_DATA,
            typename ARRAY_1D,
            typename ARRAY_1D_TO_CONST,
            typename ARRAY_2D >
  static HPCREACT_HOST_DEVICE void
  computeResidualAndJacobianAggregatePrimaryConcentrations( RealType const & temperature,
                                                            PARAMS_DATA const & params,
                                                            ARRAY_1D_TO_CONST const & targetAggregatePrimaryConcentrations,
                                                            ARRAY_1D_TO_CONST const & logPrimarySpeciesConcentration,
                                                            ARRAY_1D & residual,
                                                            ARRAY_2D & jacobian );
};



} // namespace bulkGeneric
} // namespace hpcReact

#if !defined(__INTELLISENSE__)
#include "EquilibriumReactionsAggregatePrimaryConcentration_impl.hpp"
#include "EquilibriumReactionsReactionExtents_impl.hpp"
#endif
