
#include "../DirectSystemSolve.hpp"
#include "common/pmpl.hpp"

#include <gtest/gtest.h>

using namespace hpcReact;
// TEST( testDirectSystemSolve, test3x3 )
// {
//   // **Define a Sample NxN Linear System**
//   double A_host[9] = 
//   {
//     1.0,  2.0,  3.0,
//     2.0, -1.0,  1.0,
//     3.0,  4.0,  5.0
//   };
//   double b_host[3] = { 14.0, 3.0, 24.0 }; // Right-hand side
//   double x_host[3]; // Solution

//   // **Allocate Memory on the GPU**
//   double *d_A, *d_b, *d_x;
//   cudaMalloc(&d_A, sizeof(A_host));
//   cudaMalloc(&d_b, sizeof(b_host));
//   cudaMalloc(&d_x, sizeof(x_host));

//   // **Copy Data to the GPU**
//   cudaMemcpy(d_A, A_host, sizeof(A_host), cudaMemcpyHostToDevice);
//   cudaMemcpy(d_b, b_host, sizeof(b_host), cudaMemcpyHostToDevice);

//   // **Launch Kernel (1 block, 1 thread per system)**
//   kernel_solveNxN_pivoted<double,3><<<1, 1>>>(d_A, d_b, d_x, 1);

//   // **Copy Result Back to Host**
//   cudaMemcpy(x_host, d_x, sizeof(x_host), cudaMemcpyDeviceToHost);

//   // **Print the Solution**
//   std::cout << "Solution: x = [" << x_host[0] << ", " << x_host[1] << ", " << x_host[2] << "]" << std::endl;

//   // **Free GPU Memory**
//   cudaFree(d_A);
//   cudaFree(d_b);
//   cudaFree(d_x);

// }

template< typename REAL_TYPE, int N >
struct LinearSystem
{
  REAL_TYPE A[N][N];
  REAL_TYPE b[N];
  REAL_TYPE x[N];
};

TEST( testDirectSystemSolve, test3x3 )
{
  // **Define a Sample NxN Linear System**

  LinearSystem< double, 3 > linearSystem 
  {
    { { 1.0,  2.0,  3.0 },
      { 2.0, -1.0,  1.0 },
      { 3.0,  4.0,  5.0 } 
    },
    { 14.0, 3.0, 24.0 }, // Right-hand side
    { 0.0, 0.0, 0.0 } // Solution
  };

  pmpl::genericKernelWrapper( 1, &linearSystem, [&]( auto * copyOfLinearSystem )
  {
    solveNxN_pivoted<double,3>( copyOfLinearSystem->A, copyOfLinearSystem->b, copyOfLinearSystem->x );
  } );

  EXPECT_NEAR( linearSystem.x[0], 0.0, std::numeric_limits< double >::epsilon()*100 );
  EXPECT_NEAR( linearSystem.x[1], 1.0, std::numeric_limits< double >::epsilon()*100 );
  EXPECT_NEAR( linearSystem.x[2], 4.0, std::numeric_limits< double >::epsilon()*100 );
//  std::cout << "Solution: x = [" << linearSystem.x[0] << ", " << linearSystem.x[1] << ", " << linearSystem.x[2] << "]" << std::endl;

}


#if 0 
TEST( testDirectSystemSolve, test3x3_CUDA )
{
  // **Define a Sample NxN Linear System**
  double A_host[9] = 
  {
    1.0,  2.0,  3.0,
    2.0, -1.0,  1.0,
    3.0,  4.0,  5.0
  };
  double b_host[3] = { 14.0, 3.0, 24.0 }; // Right-hand side
  double x_host[3]; // Solution

  // **Allocate Memory on the GPU**
  double *d_A, *d_b, *d_x;
  cudaMalloc(&d_A, sizeof(A_host));
  cudaMalloc(&d_b, sizeof(b_host));
  cudaMalloc(&d_x, sizeof(x_host));

  // **Copy Data to the GPU**
  cudaMemcpy(d_A, A_host, sizeof(A_host), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, b_host, sizeof(b_host), cudaMemcpyHostToDevice);

  // **Launch Kernel (1 block, 1 thread per system)**
  kernel_solveNxN_pivoted<double,3><<<1, 1>>>(d_A, d_b, d_x, 1);

  // **Copy Result Back to Host**
  cudaMemcpy(x_host, d_x, sizeof(x_host), cudaMemcpyDeviceToHost);

  // **Print the Solution**
  std::cout << "Solution: x = [" << x_host[0] << ", " << x_host[1] << ", " << x_host[2] << "]" << std::endl;

  // **Free GPU Memory**
  cudaFree(d_A);
  cudaFree(d_b);
  cudaFree(d_x);

}
#endif


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  return result;
}
