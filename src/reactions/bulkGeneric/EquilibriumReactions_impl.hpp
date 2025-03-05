#include <iostream>


namespace hpcReact
{
namespace bulkGeneric
{

constexpr bool debugPrinting = false;

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int RESIDUAL_FORM >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST,
          typename ARRAY_1D_TO_CONST2,
          typename ARRAY_2D >
HPCREACT_HOST_DEVICE inline
void
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE,
                      RESIDUAL_FORM >::computeResidualAndJacobian( RealType const & temperature,
                                                                   PARAMS_DATA const & params,
                                                                   ARRAY_1D_TO_CONST const & speciesConcentration0,
                                                                   ARRAY_1D_TO_CONST2 const & xi,
                                                                   ARRAY_1D & residual,
                                                                   ARRAY_2D & jacobian )
{

  // printf( "speciesConcentration0 = | %g %g %g %g %g |\n",
  //         speciesConcentration0[0],
  //         speciesConcentration0[1],
  //         speciesConcentration0[2],
  //         speciesConcentration0[3],
  //         speciesConcentration0[4] );
  // printf( "xi = | %g %g |\n", xi[0], xi[1] );


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
      // if( speciesConcentration[i] < 1.0e-17 )
      // {
      //   speciesConcentration[i] = 1.0e-17;
      // }
    }
  }


  // loop over reactions
  for( IndexType a=0; a<numReactions; ++a )
  {
    // get the equilibrium constant for this reaction
    RealType const Keq = params.equilibriumConstant( a );


    RealType forwardProduct = 1.0;
    RealType reverseProduct = 1.0;
    RealType dForwardProduct_dxi[numReactions] = {0.0};
    RealType dReverseProduct_dxi[numReactions] = {0.0};
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
          dForwardProduct_dxi[b] += -s_ai / speciesConcentration[i] * params.stoichiometricMatrix( b, i );
        }
      }
      else if( s_ai > 0.0 )
      {
        // reverse reaction
        reverseProduct *= pow( speciesConcentration[i], s_ai );
        // derivative of reverse product with respect to xi
        for( IndexType b=0; b<numReactions; ++b )
        {
          dReverseProduct_dxi[b] += s_ai / speciesConcentration[i] * params.stoichiometricMatrix( b, i );
        }
      }
    }
    // compute the residual for this reaction
    if constexpr( RESIDUAL_FORM == 0 )
    {
      residual[a] = forwardProduct - Keq * reverseProduct;
    }
    else if constexpr( RESIDUAL_FORM == 1 )
    {
      residual[a] = 1.0 - Keq * reverseProduct / forwardProduct;
    }
    else if constexpr( RESIDUAL_FORM == 2 )
    {
      residual[a] = log( Keq * reverseProduct / forwardProduct );
    }
    //     printf( "residual[%d] = %g - %g * %g = %g\n", a, forwardProduct, Keq, reverseProduct, residual[a] );

    // Finish the derivatives of the product terms with respect to xi
    for( IndexType b=0; b<numReactions; ++b )
    {
      // printf( "(%d) dForwardProduct_dxi[%d] = %f\n", a, b, dForwardProduct_dxi[b] );
      // printf( "(%d) dReverseProduct_dxi[%d] = %f\n", a, b, dReverseProduct_dxi[b] );
      // printf( "Keq = %f\n", Keq );

      dForwardProduct_dxi[b] *= forwardProduct;
      dReverseProduct_dxi[b] *= reverseProduct;

      if constexpr( RESIDUAL_FORM == 0 )
      {
        jacobian( a, b ) = dForwardProduct_dxi[b] - Keq * dReverseProduct_dxi[b];
      }
      else if constexpr( RESIDUAL_FORM == 1 )
      {
        jacobian( a, b ) =  -Keq * ( dReverseProduct_dxi[b] / forwardProduct - dForwardProduct_dxi[b] * reverseProduct / ( forwardProduct * forwardProduct ) );
      }
      else if constexpr( RESIDUAL_FORM == 2 )
      {
        jacobian( a, b ) = -dForwardProduct_dxi[b] / forwardProduct + dReverseProduct_dxi[b] / reverseProduct;
      }
    }
  }
}

template< typename REAL_TYPE,
          typename INT_TYPE,
          typename INDEX_TYPE,
          int RESIDUAL_FORM >
template< typename PARAMS_DATA,
          typename ARRAY_1D,
          typename ARRAY_1D_TO_CONST >
HPCREACT_HOST_DEVICE inline
void
EquilibriumReactions< REAL_TYPE,
                      INT_TYPE,
                      INDEX_TYPE,
                      RESIDUAL_FORM >::enforceEquilibrium( RealType const & temperature,
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
  for( int k=0; k<100; ++k )
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
    if( residualNorm < 1.0e-14 )
    {
      printf( " converged\n" );
      break;
    }

    if constexpr( debugPrinting )
    {
      printf( "************************************************** \n" );
      printf( "R = { " );
      for( int i=0; i<numReactions; ++i )
      {
        printf( "%18.12g", residual[i] );
        if( i < numReactions-1 )
        {
          printf( ", " );
        }
      }
      printf( " }\n" );

      printf( "J = { \n" );
      for( int i=0; i<numReactions; ++i )
      {
        printf( " { " );
        for( int j=0; j<numReactions; ++j )
        {
          printf( "%22.14g", jacobian( i, j ) );
          if( j < numReactions-1 )
          {
            printf( ", " );
          }
        }
        printf( "},\n" );
      }
      printf( " }\n" );
    }

    int numNonSymmetric = 0;
    for( int i=0; i<numReactions; ++i )
    {
      for( int j=i; j<numReactions; ++j )
      {
        if( fabs( jacobian( i, j ) - jacobian( j, i ) ) > 0.5 * ( jacobian( i, j ) + jacobian( j, i ) )*1.0e-14 )
        {
          ++numNonSymmetric;
          printf( "jacobian not symmetric: i = %d, j = %d, jacobian(i,j) = %g, jacobian(j,i) = %g\n", i, j, jacobian( i, j ), jacobian( j, i ) );
        }
      }
    }
    if( numNonSymmetric > 0 )
    {
      printf( "numNonSymmetric = %d\n", numNonSymmetric );
    }

    printf( " is J PD = %d\n", isPositiveDefinite< double, numReactions >( jacobian.data ) );

    // solve for the change in xi
    solveNxN_Cholesky< double, numReactions >( jacobian.data, residual, dxi );

    if constexpr( debugPrinting )
    {
      printf( "dxi = { " );
      for( int i=0; i<numReactions; ++i )
      {
        printf( "%18.12g", dxi[i] );
        if( i < numReactions-1 )
        {
          printf( ", " );
        }
      }
      printf( " }\n" );
    }

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
        if( cn-dc < 1.0e-30 )
        {
          REAL_TYPE const fscale = ( 1.0e-30 - cn ) / (-dc);
          if( fscale < scale )
          {
            scale = 0.99*fscale;
          }
        }
      }
    }

    if constexpr( debugPrinting )
    {
      printf( "scale = %18.12g\n", scale );
    }

    for( IndexType r=0; r<numReactions; ++r )
    {
      xi[r] -= scale * dxi[r];
    }


    if constexpr( debugPrinting )
    {
      printf( "c = { " );
      for( IndexType i=0; i<numSpecies; ++i )
      {
        REAL_TYPE c = speciesConcentration0[i];
        for( IndexType r=0; r<numReactions; ++r )
        {
          c += params.stoichiometricMatrix( r, i ) * xi[r];
        }
        printf( "%18.12g", c );
        if( i < numSpecies-1 )
        {
          printf( ", " );
        }
      }
      printf( " }\n" );
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
