#if defined(__INTELLISENSE__)
#include "EquilibriumReactions.hpp"
#endif


namespace hpcReact
{
namespace bulkGeneric
{


constexpr bool debugPrinting = false;

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline
void
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE >::computeResidualAndJacobian( REAL_TYPE const & temperature,
                                                                PARAMS_DATA const & params,
                                                                ARRAY_1D_TO_CONST const & speciesConcentration0,
                                                                ARRAY_1D_TO_CONST2 const & xi,
                                                                ARRAY_1D & residual,
                                                                ARRAY_2D & jacobian )
{

  HPCREACT_UNUSED_VAR( temperature );
  constexpr int numSpecies = PARAMS_DATA::numSpecies;
  constexpr int numReactions = PARAMS_DATA::numReactions;

  // initialize the species concentration
  RealType speciesConcentration[numSpecies];
  for( IndexType i=0; i<numSpecies; ++i )
  {
    speciesConcentration[i] = speciesConcentration0[i];
    for( IndexType r=0; r<numReactions; ++r )
    {
      speciesConcentration[i] += params.stoichiometricMatrix( r, i ) * xi[r];
    }
  }


  // loop over reactions
  for( IndexType a=0; a<numReactions; ++a )
  {
    // get the equilibrium constant for this reaction
    RealType const Keq = params.equilibriumConstant( a );


    RealType forwardProduct = 1.0;
    RealType reverseProduct = 1.0;

    // these actually only hold the derivatives of the single concentration term for the product...not the full product.
    // it will have to be multiplied by the product itself to get the derivative of the product.
    RealType dForwardProduct_dxi_divProduct[numReactions] = {0.0};
    RealType dReverseProduct_dxi_divProduct[numReactions] = {0.0};
    // loop over species
    for( IndexType i=0; i<numSpecies; ++i )
    {
      RealType const s_ai = params.stoichiometricMatrix( a, i );
      if( s_ai < 0.0 )
      {
        // forward reaction
        forwardProduct *= pow( speciesConcentration[i], -s_ai );
        // derivative of forward product with respect to xi
        for( IndexType b=0; b<numReactions; ++b )
        {
          dForwardProduct_dxi_divProduct[b] += -s_ai / speciesConcentration[i] * params.stoichiometricMatrix( b, i );
        }
      }
      else if( s_ai > 0.0 )
      {
        // reverse reaction
        reverseProduct *= pow( speciesConcentration[i], s_ai );
        // derivative of reverse product with respect to xi
        for( IndexType b=0; b<numReactions; ++b )
        {
          dReverseProduct_dxi_divProduct[b] += s_ai / speciesConcentration[i] * params.stoichiometricMatrix( b, i );
        }
      }
    }
    // compute the residual for this reaction
    residual[a] = log( reverseProduct / ( forwardProduct * Keq ) );

    // compute the jacobian
    for( IndexType b=0; b<numReactions; ++b )
    {
      jacobian( a, b ) = -dForwardProduct_dxi_divProduct[b] + dReverseProduct_dxi_divProduct[b];
    }
  }
}



template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST >
HPCREACT_HOST_DEVICE inline
void
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE >::enforceEquilibrium_Extents( REAL_TYPE const & temperature,
                                                        PARAMS_DATA const & params,
                                                        ARRAY_1D_TO_CONST const & speciesConcentration0,
                                                        ARRAY_1D & speciesConcentration )
{
  HPCREACT_UNUSED_VAR( temperature );
  constexpr int numSpecies = PARAMS_DATA::numSpecies;
  constexpr int numReactions = PARAMS_DATA::numReactions;
  double residual[numReactions] = { 0.0 };
  double xi[numReactions] = { 0.0 };
  double dxi[numReactions] = { 0.0 };
  CArrayWrapper< double, numReactions, numReactions > jacobian;

  REAL_TYPE residualNorm = 0.0;
  for( int k=0; k<30; ++k )
  {
    computeResidualAndJacobian( temperature,
                                params,
                                speciesConcentration0,
                                xi,
                                residual,
                                jacobian );

    residualNorm = 0.0;
    for( int j = 0; j < numReactions; ++j )
    {
      residualNorm += residual[j] * residual[j];
    }
    residualNorm = sqrt( residualNorm );
    printf( "iter, residualNorm = %2d, %16.10g \n", k, residualNorm );
    if( residualNorm < 1.0e-12 )
    {
      printf( " converged\n" );
      break;
    }

    // solve for the change in xi
    for( int r=0; r<numReactions; ++r )
    {
      dxi[r] = 0.0;
      residual[r] = -residual[r];
    }

    solveNxN_Cholesky< double, numReactions >( jacobian.data, residual, dxi );


    // scaling
    REAL_TYPE scale = 1.0;
    for( IndexType i=0; i<numSpecies; ++i )
    {
      REAL_TYPE cn = speciesConcentration0[i];
      REAL_TYPE dc = 0.0;
      for( IndexType r=0; r<numReactions; ++r )
      {
        dc += params.stoichiometricMatrix( r, i ) * dxi[r];
        cn += params.stoichiometricMatrix( r, i ) * xi[r];
      }
      if( cn+dc < 1.0e-30 )
      {
        REAL_TYPE const fscale = ( 1.0e-30 - cn ) / (dc);
        if( fscale < scale )
        {
          scale = 0.9*fscale;
          //printf( "i, cn, dc, scale = %d, %18.12g %18.12g %18.12g\n", i, cn, dc, scale );
        }
      }
    }

    for( IndexType r=0; r<numReactions; ++r )
    {
      xi[r] += scale * dxi[r];
    }

  }

  for( IndexType i=0; i<numSpecies; ++i )
  {
    speciesConcentration[i] = speciesConcentration0[i];
    for( IndexType r=0; r<numReactions; ++r )
    {
      speciesConcentration[i] += params.stoichiometricMatrix( r, i ) * xi[r];
    }
  }
}


} // namespace bulkGeneric
} // namespace hpcReact