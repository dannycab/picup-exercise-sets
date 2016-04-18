#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main () {
	
  /* Declare variables */
	int tsteps,j;
  double k,m,c,vi,xi,dt;
  double *t,*x,*v,*T,*U,*E;
  FILE *v_file,*x_file,*E_file;

  k = 10;
  m = 0.01;
  c = 0.02;
  vi = 0;
  xi = 1;
  dt = 0.001;
  tsteps = 3000;

  /* Allocate memory and initialize contents to zero */
  t=(double*)calloc(tsteps+2, sizeof(double));
  x=(double*)calloc(tsteps+2, sizeof(double));
  v=(double*)calloc(tsteps+2, sizeof(double));
  T=(double*)calloc(tsteps+2, sizeof(double));
  U=(double*)calloc(tsteps+2, sizeof(double));
  E=(double*)calloc(tsteps+2, sizeof(double));

  x[1] = xi;
  v[1] = vi;
  t[1] = 0;

  /* Open the files used to store the results */
	v_file=fopen("vFile.dat","w");
	x_file=fopen("xFile.dat","w");
	E_file=fopen("Efile.dat","w");

	for(j=1; j<=tsteps; j++) {
		
    t[j]=((double)j)*dt;
    v[j+1] = v[j] - (k/m) * x[j] * dt - (c/m) * v[j] * dt;
    x[j+1] = x[j] + v[j+1] * dt;

    T[j+1] = 0.5*pow(v[j+1],2);
    U[j+1] = 0.5*k*pow(x[j+1],2);
    E[j+1] = T[j+1] + U[j+1];

		fprintf(v_file,"%lf\t%10.14G\n",t[j],v[j+1]);
		fprintf(x_file,"%lf\t%10.14G\n",t[j],x[j+1]);
		fprintf(E_file,"%lf\t%10.14G\t%10.14G\t%10.14G\n",t[j],T[j+1],U[j+1],E[j+1]);
	}
	
  /* Deallocate memory */
  free(t);
  free(v);
  free(x);
  free(T);
  free(U);
  free(E);

  /* Close the files */
	fclose(v_file);
	fclose(x_file);
	fclose(E_file);

  return 0;
}