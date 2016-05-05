#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Relaxation method (Explicit)

int main () {

    // Declare variables, mesh, and output file

    int n,m,k,imax;
    double L,W,Dx,Dy,R,beta,denom,epsilon;
    double *x,*y;
    double **X, **Y, **V,**VN,**err;
    FILE *V_file;

    // Initialize variables needed for mesh
    L = 5.0;
    W = 5.0;
    n = 10;
    m = 10;
    Dx = L/n;
    Dy = W/m;

    // Allocate memory for spatial vectors
    x=(double*)calloc(n, sizeof(double));
    y=(double*)calloc(m, sizeof(double));

    // Allocate memory for matrices
    X=(double**) malloc(n, sizeof(double));
    Y=(double**) malloc(n, sizeof(double));
    V=(double**) malloc(n, sizeof(double));
    VN=(double**) malloc(n, sizeof(double));
    err=(double**) malloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {

        X[i] = (double*) malloc(m, sizeof(double));
        Y[i] = (double*) malloc(m, sizeof(double));
        V[i] = (double*) malloc(m, sizeof(double));
        VN[i] = (double*) malloc(m, sizeof(double));
        err[i] = (double*) malloc(m, sizeof(double));

    }


    // Open file for output
    V_file = fopen("VFile.dat","w");

    R = 5.0;

    //Is there an equivalent ones function in C?
    //V = R*np.ones((n+1,m+1));

    for(i=1; j<=n; j++) { //This might be wrong because of 0 vs 1 counting in C vs Python

        V[i][m] = 10;

    }

    VN = V; // New iterated solution; deals with reference issue
    err = VN - V; // Error between old and new iterated solution

    beta = Dx/Dy; // parameter from finite difference algorithm
    denom = 2*(1+beta**2); // just to make the denominator easy

    // Iterative procedure (tolerances and counters) 

    epsilon = 1e-5;
    imax = 5000;
    k = 1;

    while k < imax:
        
        for i in np.arange(1,n):
            
            for j in np.arange(1,m):
                
                VN[i,j]=(V[i-1,j]+V[i+1,j]+beta**2*(V[i,j-1]+V[i,j+1]))/denom;
                err[i,j] = VN[i,j]-V[i,j];
        
        V = VN + np.zeros((n+1,m+1)); # Deals with reference issue
        k += 1;
        
        errmax = np.max(err); 
               
        if errmax < epsilon:
            
            print("Convergence after", k, "iterations.");
            break;
            
        if k == imax:
            
            print("Did not converge after", k, "iterations.");
            print("Max error on mesh: ", errmax);



    // Free memory solution written to disk

    for (i = 0; i < n; i++){

        free(X[i]);
        free(Y[i]);
        free(V[i]);
        free(VN[i]);
        free(err[i]);

    }

    free(X);
    free(Y);
    free(V);
    free(VN);
    free(err);


    fclose(V_file);

}


