
#pragma once

#include "common/macros.hpp"

#include <cmath>

namespace hpcReact
{
namespace bulkGeneric
{


template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE
          >
class ReactionsBase
{
public:
  using RealType = REAL_TYPE;
  using RealDataArrayView1d = REAL_DATA_ARRAY_1D_VIEW_TYPE;
  using RealConstDataArrayView1d = REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE;
  using IntType = INT_TYPE;
  using IndexType = INDEX_TYPE;


  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE void
  computeLog10ActCoefBDotModel( REAL_TYPE const temperature,
                                REAL_TYPE const ionicStrength,
                                PARAMS_DATA const & params,
                                RealDataArrayView1d & log10PrimaryActCoeff,
                                RealDataArrayView1d & dLog10PrimaryActCoeff_dIonicStrength,
                                RealDataArrayView1d & log10SecActCoeff,
                                RealDataArrayView1d & dLog10SecActCoeff_dIonicStrength );

  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE void
  computeIonicStrength( PARAMS_DATA const & params,
                        RealConstDataArrayView1d const & primarySpeciesConcentration,
                        RealConstDataArrayView1d const & secondarySpeciesConcentration,
                        REAL_TYPE & ionicStrength );
};

} // namespace bulkGeneric
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
