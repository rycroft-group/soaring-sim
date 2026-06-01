#include "common.hh"
#include "lp_solve.hh"
#include <cstdio>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <sys/stat.h>

// Test matrix update code
int main(){

	const int n=10; // size of matrix
	const int vm=2; // index of deleted column/row

	double* A=new double[n*n];double* A_inv=new double[n*n];
    double* R=new double[n];double* S=new double[n];double* V=new double[n];
    gsl_rng *rng;
    rng=gsl_rng_alloc(gsl_rng_taus2);
    // Initialize A
    for(int j=0;j<n;j++){
    	for(int i=0;i<=j;i++){
    		double val=gsl_rng_uniform(rng);
    		A[j*n+i]=val;
    		A_inv[j*n+i]=val;
		printf("j %d i %d A %g\n",j,i,A[j*n+i]);

    	}
    }
	matrix_inv(A_inv,n);
    // Compute R=A^-1 e
    // R is the vth row/column of A
    for(int i=0;i<n;i++){
        if(i<vm){
            R[i]=A_inv[i+vm*n];
        }
        else{
            R[i]=A_inv[vm+i*n];
        }
    }
    // Compute S=A^-1 V
    // Init V, the changed values
    for(int i=0;i<n;i++){
//	    printf("old M %g %g %g\n",old_M.x,M[vm].y,M[i].z);
        V[i]=gsl_rng_uniform(rng)+1.;
    }
    V[vm]=0.;

    for(int i=0;i<n;i++){
       S[i]=0.;
       for(int j=0;j<n;j++){
            if(j<i){
                S[i]+=A_inv[j+i*n]*V[j];
            }
            else{
                S[i]+=A_inv[i+j*n]*V[j];
            }
        }
    }
    // Copute alpha,beta,gamma and delta
    double alpha=R[vm];
    double beta=1.+S[vm];
    double gamma=0.;
    for(int i=0;i<n;i++){
        gamma+=V[i]*S[i];
    }
    double sigma=1./(beta*beta-alpha*gamma);
    // update A^{-1}
     // Create covariance matrix
    for(int j=0;j<n;j++){
        for(int i=0;i<=j;i++){
            A_inv[j*n+i]+=gamma*sigma*R[i]*R[j]-beta*sigma*(R[i]*S[j]+S[i]*R[j])+alpha*sigma*S[i]*S[j];
            // A[j*n+i]=0.1*i;
        }
    }

    for(int j=0;j<n;j++){
    	if(j<vm){
    		A[vm*n+j]+=V[j];
    	}
    	else{
    		A[j*n+vm]+=V[j];
    	}

    }
    for(int j=0;j<n;j++)
           for(int i=0;i<=j;i++){
                   printf("j %d i %d update A %g\n",j,i,A[j*n+i]);
           }
    matrix_inv(A,n);

    for(int j=0;j<n;j++)
	   for(int i=0;i<=j;i++){
		   printf("j %d i %d inv A %g Solution %g\n",j,i,A_inv[j*n+i],A[j*n+i]);
	   }
    delete [] A;
    delete [] V;
    delete [] S;
    delete [] R;
    delete [] A_inv;

}
