#include "macros.hpp"

#include<cmath>

namespace hpcReact
{
template< typename REAL_TYPE, int N >
HPCREACT_HOST_DEVICE void solveNxN_pivoted(REAL_TYPE A[N][N], REAL_TYPE b[N], REAL_TYPE x[N]) 
{
    int pivot[N] = {0, 1, 2};  // Row index tracker
    for (int i = 0; i < N; i++) 
    {
        pivot[i] = i;
    }

    // **Step 1: Forward Elimination with Pivoting**
    for (int k = 0; k < N-1; k++) 
    {
        // **Find Pivot Row**
        int max_row = k;
        REAL_TYPE max_val = abs(A[pivot[k]][k]);
        for (int i = k + 1; i < N; i++) 
        {
            if (fabs(A[pivot[i]][k]) > max_val) 
            {
                max_val = fabs(A[pivot[i]][k]);
                max_row = i;
            }
        }

        // **Swap Rows in Pivot Array**
        if (max_row != k) 
        {
            int temp = pivot[k];
            pivot[k] = pivot[max_row];
            pivot[max_row] = temp;
        }

        // **Gaussian Elimination**
        for (int i = k + 1; i < N; i++) 
        {
            REAL_TYPE factor = A[pivot[i]][k] / A[pivot[k]][k];
            for (int j = k; j < N; j++) 
            {
                A[pivot[i]][j] -= factor * A[pivot[k]][j];
            }
            b[pivot[i]] -= factor * b[pivot[k]];
        }
    }

    // **Step 2: Back-Substitution**
    for (int i = N - 1; i >= 0; --i) {
        x[i] = b[pivot[i]];
        for (int j = i + 1; j < N; j++) {
            x[i] -= A[pivot[i]][j] * x[j];
        }
        x[i] /= A[pivot[i]][i]; // Normalize
}

}

} // namespace hpcReact