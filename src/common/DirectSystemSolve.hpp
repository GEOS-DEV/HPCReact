#include "macros.hpp"

#include<cmath>
template< typename REAL_TYPE, int N >
HPCREACT_HOST_DEVICE void solveNxN_pivoted(REAL_TYPE A[N][N], REAL_TYPE b[N], REAL_TYPE x[N]) 
{
    int pivot[3] = {0, 1, 2};  // Row index tracker

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
    x[2] = b[pivot[2]] / A[pivot[2]][2];
    x[1] = (b[pivot[1]] - A[pivot[1]][2] * x[2]) / A[pivot[1]][1];
    x[0] = (b[pivot[0]] - A[pivot[0]][1] * x[1] - A[pivot[0]][2] * x[2]) / A[pivot[0]][0];
}

