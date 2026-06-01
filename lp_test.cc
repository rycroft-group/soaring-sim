#include <cmath>
#include <cstdio>
#include "lp_solve.hh"

int main(){

	// test symmetric matrix-vector mul
	const int m=4,n=4;
	double* A=new double[m*n];
	puts("A before mul");
	for(int j=0;j<n;j++){
		for(int i=0;i<j;i++){
			A[j*n+i]=sin(i)*j;
			printf("i %d, j %d, A[ij] %f\n",i,j,A[j*n+i]);
		}
		A[j*n+j]=1.;
	}// correct

	double B[4*2]={1.,2.,3.,4.,5,6,7,8};
	double* C=new double[4*2];
	matrix_mul(A,B,C,m,n,2,true);
	puts("A after mul");
	for(int j=0;j<n;j++){
                for(int i=0;i<j;i++){
                        // A[j*n+i]=1.*i*j;
                        printf("i %d, j %d, A[ij] %f\n",i,j,A[j*n+i]);
                }
                // A[j*n+j]=1.;
        }//
	for(int i=0;i<n*2;i++)
		printf("i %d,C[i] %f\n",i,C[i]);

	// test symmetric matrix inverse
	matrix_inv(A,n);
	for(int j=0;j<n;j++){
		for(int i=0;i<j;i++){
			printf("i %d, j %d, A[ij] %f\n",i,j,A[j*n+i]);
		}
		printf("i %d, j %d, A[ij] %f\n",j,j,A[j*n+j]);
	} // correct

	// test matrix update
	const int mm=m-1,nn=n-1;
	double P[m*n];
	for(int j=0;j<nn;j++){
		for(int i=0;i<j;i++){
			puts("?");
			P[j*nn+i]=sin(i)*j;
			printf("i %d, j %d, P[ij] %f\n",i,j,P[j*nn+i]);
		}
		P[j*nn+j]=1.;
	}
	puts("inverse");
	matrix_inv(P,nn);
	for(int j=0;j<nn;j++){
		for(int i=0;i<j;i++){
			printf("i %d, j %d, P[ij] %f\n",i,j,P[j*nn+i]);
		}
		printf("i %d, j %d, P[ij] %f\n",j,j,P[j*nn+j]);

	}

	B[0]=0.,B[1]=sin(1)*3.,B[2]=sin(2.)*3.;
	double D=1.;
	puts("update");
	matrix_update(P,B,D,n);
	for(int j=0;j<n;j++){
		for(int i=0;i<j;i++){
			printf("i %d, j %d, P[ij] %f\n",i,j,P[j*n+i]);
		}
		printf("i %d, j %d, P[ij] %f\n",j,j,P[j*n+j]);
	}

	matrix_mul(P,C,B,n,n,2,true);
	for(int i=0;i<n*2;i++) printf("i %d B[i] %f\n",i,B[i]);//correct

	// test B=A^{-1}C is equivalent to solve(AB=C)
	puts("test solve");
	for(int j=0;j<n;j++){
		for(int i=0;i<j;i++){
			A[j*n+i]=sin(i)*j;
			printf("i %d, j %d, A[ij] %f\n",i,j,A[j*n+i]);
		}
		A[j*n+j]=1.;
		printf("i %d, j %d, A[ij] %f\n",j,j,A[j*n+j]);
	}
	for(int i=0;i<2*n;i++){
		B[i]=1.*i;
		printf("i %d B[i] %f\n",i,B[i]);
	}
	matrix_mul(A,B,C,n,n,2,true);
	for(int i=0;i<2*n;i++) printf("i %d,C[i] %f\n",i,C[i]);
	int nrhs=2;
	solve_sym_matrix(n,A,C,nrhs);
	for(int i=0;i<2*n;i++) printf("i %d,B[i] %f\n",i,C[i]);

	matrix_mul(A,B,C,n,n,2,true);
	matrix_inv(A,n);
	matrix_mul(A,C,B,n,n,2,true);
	for(int i=0;i<2*n;i++) printf("i %d,B[i] %f\n",i,B[i]);

	delete [] C;
	delete [] A;

}
