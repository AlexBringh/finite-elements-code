//
//	Author: Alexander B. Ringheim
//	Project: Finite Element Method solver for Plasticity analysis of materials and structures.
//	Date of creation: 14.10.2024
//
#include "gauss_elimination.h"
#include "matrix_utils.h"

int gauss_elimination (double *A, int m, double *B)
{
    /*
        Performs Gauss Elimination on a augmented matrix.
        The matrix must be m-rows and n = m + 1 columns.
    */
    int forw = forward_elimination(A, m);
    int back = backwards_substitution(A, m, B);

    return 0;
}

int forward_elimination (double *A, int m)
{
    int n = m + 1;
    printf("Forwards elimination. \nInput matrix: \n");
    print_matrix(A, m, m+1);
    for (int i = 0; i < m; i++)
    {
        // Perform partial pivoting
        int max = i;
        for (int k = i + 1; k < m; k++) // i + 1 references to the columns.
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
            for (int j = 0; j <= m; j++)
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
        for (int j = 0; j <= m; j++)
        {
            *(A + i * n + j) /= p;
        }

        // Eliminate the entries below the diagonal.
        for (int j = i + 1; j < m; j++)
        {
            double f = *(A + j * n + i); 
            for (int k = 0; k <= m; k++)
            {
                *(A + j * n + k) -= f * *(A + i * n + k);
            }
        }


        printf("Forwards elimination, step: %d \n", i);
        print_matrix(A, m, n);
    }

    return 0;
}


int backwards_substitution (double *A, int m, double *B)
{
    int n = m + 1;
    for (int i = m - 1; i >= 0; --i)
    {
        *(B + i) = *(A + i * n + m);
        for (int j = i + 1; j < m; j++)
        {
            *(B + i) -= *(A + i * n + j) * *(B + j);
        }
    }

    printf("Resulting matrix: \n");
    print_matrix(B, m, 1);

    return 0;
}