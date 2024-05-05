#include <stdio.h>
#include <stdlib.h>

void initialize(int n, double A[n][n], double a[n], double b[n], double c[n], double d[n]);

int main(){

    int i, j, k, n;

    printf("Enter size of Matrix (n): ");
    scanf("%d", &n);
    double A[n][n+1];
    double a[n], b[n], c[n], d[n];
    double x[n], sum = 0.0;

    // Initialize
    initialize(n, A, a, b, c, d);

    // Print Matrix A
    printf("Matrix A:\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j <= n; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }

    // generation of upper triangular matrix
    for (j = 0; j < n - 1; j++) {
        for (int i = j + 1; i < n; i++) {
            double w = A[i][j] / A[j][j];
            for (int k = j; k <= n; k++) {
                A[i][k] = A[i][k] - w * A[j][k];
            }
        }
    }

    // Print Upper Triangular Matrix
    printf("\nUpper Triangular Matrix A:\n");
    for (i = 0; i < n; i++) {
        for (j = 0; j <= n; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }

    /* backward substitution*/
    x[n-1] = A[n-1][n] / A[n-1][n-1];
    for (i = n - 2; i >= 0; i--)
    {
        sum = 0;
        for (j = i + 1; j < n; j++)
        {
            sum = sum + A[i][j] * x[j];
        }
        x[i] = (A[i][n-1] - sum) / A[i][i];
    }

    printf("\nThe solution is: \n");
    for (i = 0; i < n; i++)
    {
        printf("\nx%d=%2f\t", i, x[i]);
    }
    printf("\n");
}


void initialize(int n, double A[n][n+1], double a[n], double b[n], double c[n], double d[n]) {

    int i, j;
    for (i=0; i<n; i++){
        a[i] = 4;
        d[i] = 1;
    }

    for (i=1; i<n; i++){
        b[i] = 1;
    }

    for (i=0; i<n-1; i++){
        c[i] = 2;
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j <= n; j++) {
            if (i == j) {
                A[i][j] = a[i];
            } else if (j == i - 1) {
                A[i][j] = b[i];
            } else if (j == i + 1) {
                A[i][j] = c[i];
            } else {
                A[i][j] = 0;
            }
        }
        A[i][n] = d[i];
    }
}

