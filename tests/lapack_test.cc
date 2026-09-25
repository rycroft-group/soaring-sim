#include <cstdio>
#include <cstdlib>

// Tell the compiler about the existence of the required LAPACK functions
extern "C" {
    void dsytrf_(char* uplo,int *n,double *a,int *lda,int *ipivot,double *work,
            int *lwork,int *info);
    void dsytri_(char* uplo,int *n,double *A,int *lda,int *ipivot,double *work,int *info);
}

int main() {
    double A[9]={1,0,0,0.2,2,0,-1.,0,3};
    double b[3]={3,7,2};

    // Print the matrix and the right hand side
    printf("A=[%6g %6g %6g]\n  [%6g %6g %6g]\n  [%6g %6g %6g]\n\n",*A,A[3],A[6],A[1],A[4],A[7],A[2],A[5],A[8]);
    //printf("b=[%6g ]\n  [%6g ]\n  [%6g ]\n\n",*b,b[1],b[2]);

    // Call routine to solve the matrix
    char uplo='u';
    int info,n=2,m=3;
    int lwork=1024;
    int ipiv[3];
    double work[1024];

    // Perform LU decomposition
    dsytrf_(&uplo,&n,A,&m,ipiv,work,&lwork,&info);

    // Perform inversion
    dsytri_(&uplo,&n,A,&m,ipiv,work,&info);

    printf("A=[%6g %6g %6g]\n  [%6g %6g %6g]\n  [%6g %6g %6g]\n\n",*A,A[3],A[6],A[1],A[4],A[7],A[2],A[5],A[8]);
    // Print the solution
    //printf("x=[%6g ]\n  [%6g ]\n",*b,b[1],b[2]);
}
