#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Relaxation method (Explicit)

int main () {

    // Declare variables and output file

    int n,m,k,imax;
    double L,W,Dx,Dy,R,beta,denom,epsilon;
    double *x,*y;
    FILE *V_file;

    // Initialize variables needed for mesh
    L = 5.0;
    W = 5.0;
    n = 50;
    m = 50;
    Dx = L/n;
    Dy = W/m;

    // Declare mesh
    // X,Y = np.meshgrid(x,y); Need to finish this
    double X[n][m],Y[n][m],V[m][n],VN[m][n],err[m][n];

    // Allocate memory for spatial vectors
    x=(double*)calloc(n, sizeof(double));
    y=(double*)calloc(m, sizeof(double));

    // How do I allocate memory for a matrix?

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
        
        #print(err)
        
        if errmax < epsilon:
            
            print("Convergence after", k, "iterations.");
            break;
            
        if k == imax:
            
            print("Did not converge after", k, "iterations.");
            print("Max error on mesh: ", errmax);


// plt.figure()
// CS = plt.contour(X, Y, V.transpose())
// plt.title('Potential contours')
// plt.xlim((0,L));
// plt.ylim((0,W));
// CB = plt.colorbar(CS, shrink=0.8, extend='both');

// # 3D Plot
// fig = plt.figure()
// ax = fig.gca(projection='3d')
// surf = ax.plot_surface(X, Y, V, rstride=1, cstride=1, cmap=cm.coolwarm,
//                        linewidth=0, antialiased=False)
// fig.colorbar(surf, shrink=0.5, aspect=5)

// ax.set_xlabel('x')
// ax.set_ylabel('y')
// ax.set_zlabel('Potential (V)')
// plt.show()



