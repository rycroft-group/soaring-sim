#include <cstdio>
#include <cstdlib>

#include <cstdio>
#include <cstdlib>

// Tell the compiler about the existence of the required LAPACK functions
extern "C" {
    int dgetrs_(char *trans_,int *n,int *nrhs,double *a,int *lda,
            int *ipiv,double *b,int *ldb,int *info);
    int dgetrf_(int *m,int *n,double *a,int *lda,int *ipiv,int *info);
}

// Solves the matrix system Ax=b, returning the answer in the b array
void solve_matrix(int n,double *A,double *b) {

    // Create the temporary memory that LAPACK needs
    int info,nrhs=1,*ipiv=new int[n];
    char trans='N';

    // Perform the LU decomposition
    dgetrf_(&n,&n,A,&n,ipiv,&info);
    if(info!=0) {
        fputs("LAPACK LU routine failed\n",stderr);
        exit(1);
    }

    // Use the LU decomposition to solve the system
    dgetrs_(&trans,&n,&nrhs,A,&n,ipiv,b,&n,&info);
    if(info!=0) {
        fputs("LAPACK solve routine failed\n",stderr);
        exit(1);
    }

    // Remove temporary memory
    delete [] ipiv;
}

int main() {

    // Create test 3x3 matrix A and right hand side vector b
    double A[9]={2,0.3,0.1,0.2,2,0,-1.,0,3};
    double b[3]={3,7,2};

    // Print the matrix and the right hand side
    printf("A=[%7.4g %7.4g %7.4g]\n  [%7.4g %7.4g %7.4g]\n  [%7.4g %7.4g %7.4g]\n\n",*A,A[3],A[6],A[1],A[4],A[7],A[2],A[5],A[8]);
    printf("b=[%7.4g ]\n  [%7.4g ]\n  [%7.4g ]\n\n",*b,b[1],b[2]);

    // Call routine to solve the matrix
    solve_matrix(2,A,b);

    // Output the modified state of A, which encodes the L and U matrices
    printf("A (modified) =[%7.4g %7.4g %7.4g]\n              "
           "[%7.4g %7.4g %7.4g]\n              [%7.4g %7.4g %7.4g]\n\n",*A,A[3],A[6],A[1],A[4],A[7],A[2],A[5],A[8]);

    // Print the solution x
    printf("x=[%7.4g ]\n  [%7.4g ]\n  [%7.4g ]\n",*b,b[1],b[2]);
}
