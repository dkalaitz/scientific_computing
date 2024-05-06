#include<stdio.h>
#include<math.h>
#include <stdlib.h>
#include<omp.h>

// Function to calculate the maximum norm of the difference between two vectors
double maxNorm(double x1[],double x2[],int n){
	double maxVal=fabs(x1[0]-x2[0]),tmp;
	int i;
	for(i=1;i<n;i++){	
		tmp=fabs(x1[i]-x2[i]);
		if(tmp>maxVal)
			maxVal=tmp;
	}
	return maxVal;
}

int main()
{
    int i,j,k,iters=0,maxIters=100,n, numThreads;
    double omega;

    // Input size of Matrix (n)
    printf("Enter size of Matrix (n): ");
    scanf("%d", &n);

    // Allocate memory for matrices and vectors
    double **A, *b, *x_old, *x_new;
    double sum=0.0,eps=0.0000001,error;

    b = (double *)malloc(n * sizeof(double));
    x_old = (double *)malloc(n * sizeof(double));
    x_new = (double *)malloc(n * sizeof(double));
    A = (double **)malloc(n * sizeof(double *));

    if (!A || !b || !x_old || !x_new) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    for (i = 0; i < n; i++) {
        A[i] = (double *)malloc(n * sizeof(double));
        if (A[i] == NULL) {
            printf("Memory allocation failed.\n");
            return 1;
        }
    }

    // Initialize Matrix A, vector x and b
    for (i=0; i<n; i++){
        x_old[i] = 0;
        b[i] = 1;
        for (j=0; j<n; j++){
            if (i==j){
                A[i][j] = 4;
            } else if (j == i-1 || j == i+1){
                A[i][j] = -1;
            } else {
                A[i][j] = 0;
            }
        }
    }

    // Print Matrix A if small
    if (n < 7){
        printf("Matrix A:\n");
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                printf("%f ", A[i][j]);
            }
            printf("\n");
        }
    }
    
    // Input number of threads and relaxation parameter (omega)
    printf("\nEnter number of threads: ");
    scanf("%d", &numThreads);
    printf("\nEnter relaxation parameter (omega): ");
    scanf("%lf", &omega);

    double start_time, end_time;
    start_time = omp_get_wtime();

    #pragma omp parallel num_threads(numThreads) shared(x_old, x_new, A, b) private(i, j, sum)
    {
        do {
            #pragma omp for 
                for(i=0;i<n;i++){
                    sum=0.0;
                    for(j=0;j<n;j++){
                        if(i!=j)
                            sum += A[i][j] * x_new[j];
                    }
                    x_old[i] = x_new[i];
                    x_new[i] = (1 - omega) * x_old[i] + (omega / A[i][i]) * (b[i] - sum);
                }
            #pragma omp single 
            {
                iters++;
                error=maxNorm(x_old,x_new,n);
                for(i=0;i<n;i++)
                    x_old[i]=x_new[i];
            }
        } while(error > eps && iters < maxIters);
    }
   	
    end_time = omp_get_wtime();
   	double cpu_time_used = (double) (end_time - start_time);
   	
   	// Print the solution vector
	printf("\nThe solution is: \n");
    for(i=0; i<n; i++)    
        printf("\nx%d=%.10f\t",i,x_new[i]);  

    // Print number of iterations and CPU time used
    printf("\n\nRequired number of iterations for convergence: %d\n",iters);
    printf("\nCPU time used for SOR: %f seconds\n", cpu_time_used);
    
    // Free dynamically allocated memory
    for (int i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    free(x_old);
    free(x_new);
}