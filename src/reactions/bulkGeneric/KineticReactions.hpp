#pragma once

#include "common/macros.hpp"

namespace hpcReact
{
namespace bulkGeneric
{

template< typename REAL_TYPE,
          typename REAL_DATA_ARRAY_1D_VIEW_TYPE,
          typename REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
class KineticReactions
{
public:

  using RealType = REAL_TYPE;
  using RealDataArrayView1d = REAL_DATA_ARRAY_1D_VIEW_TYPE;
  using RealConstDataArrayView1d = REAL_CONST_DATA_ARRAY_1D_VIEW_TYPE;
  using IntType = INT_TYPE;
  using IndexType = INDEX_TYPE;

  template< typename PARAMS_DATA >
  static HPCREACT_HOST_DEVICE inline void
  computeReactionRates( RealType const & temperature,
                        PARAMS_DATA const & params,
                        RealConstDataArrayView1d const & primarySpeciesConcentration,
                        RealDataArrayView1d & reactionRates );


};

} // namespace bulkGeneric
} // namespace hpcReact

#include "common/macrosCleanup.hpp"
