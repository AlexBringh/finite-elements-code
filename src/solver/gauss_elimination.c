//
//	Author: Alexander B. Ringheim
//	Project: Finite Element Method solver for Plasticity analysis of materials and structures.
//	Date of creation: 14.10.2024
//
#include "gauss_elimination.h"
#include "print_matrix.h"

int gauss_elimination (double *A, int n, double *B)
{
    int forw = forward_elimination(A, n);
    int back = backwards_substitution(A, n, B);

    return 0;
}

int forward_elimination (double *A, int n)
{
    printf("Forwards elimination. \nInput matrix: \n");
    print_matrix(A, n, n+1);
    for (int i = 0; i < n; i++)
    {
        // Perform partial pivoting
        int max = i;
        for (int k = i + 1; k < n; k++)
        {
            double a = *(A + k * n + i);
            double b = *(A + max * n + i);
            /* If the item in the first row at the i-th column is larger than the same column item for the row below, swap the rows.
               This will loop through for all the values, however moving along the diagonal. 
               The loop will continue till the largest possible value exists on the diagonal for each variable (however, favoring the first variables)
            */
            if (fabs(a) > fabs(b))
            {
                max = k;
            }
        }

        // Swap the row with the largest value with the current row (the pivot row)
        if (max != i)
        {
            for (int j = 0; j <= n; j++)
            {
                double t = *(A + i * n + j);
                *(A + i * n + j) = *(A + max * n + j);
                *(A + max * n + j) = t;
            }
        }

        // Make each diagonal element 1.
        double p = *(A + i * n + i);
        if (p == 0)
        {
            fprintf(stderr, "Matrix is singular or nearly singular! \n");
            return -1;
        }
        for (int j = 0; j <= n; j++)
        {
            *(A + i * n + j) /= p;
        }

        // Eliminate the entries below the diagonal.
        for (int j = i + 1; j < n; j++)
        {
            double f = *(A + j * n + i); 
            for (int k = 0; k <= n; k++)
            {
                *(A + j * n + k) -= f * *(A + i * n + k);
            }
        }


        printf("Forwards elimination, step: %d \n", i);
        print_matrix(A, n, n+1);
    }

    return 0;
}


int backwards_substitution (double *A, int n, double *B)
{
    for (int i = n - 1; i >= 0; --i)
    {
        *(B + i) = *(A + i * n + n);
        for (int j = i + 1; j < n; j++)
        {
            *(B + i) -= *(A + i * n + j) * *(B + j);
        }
    }

    printf("Resulting matrix: \n");
    print_matrix(B, n, 1);

    return 0;
}