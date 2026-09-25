#include "lp_solve.hh"
#include <cstdio>
#include <cstdlib>
#include <cstring>

// Tell the compiler about the existence of the required LAPACK functions
extern "C" {
    void dsysv_(char* uplo,int *n,int *nrhs,double *a,int *lda,int *ipivot,
            double *b,int *ldb,double *work,int *lwork,int *info);
    void dsytrf_(char* uplo,int *n,double *a,int *lda,int *ipivot,double *work,
            int *lwork,int *info);
    void dsytrs_(char* uplo,int *n,int *nrhs,double *a,int *lda,int *ipivot,
            double *b,int *ldb,int *info);
    void dsytri_(char* uplo,int *n,double *A,int *lda,int *ipivot,double *work,int *info);
    void dsymv_(char* uplo,int *n, double *alpha,double *A,int *lda,double *B,int *incx,
            double *beta,double *Y,int *incy);
    void dsymm_(char* side,char* uplo,int *m,int *n,double *alpha,double* A,int *lda,
            double *B,int *ldb,double *beta,double *C,int *ldc);
    void dgemv_(char* trans,int *m,int *n,double *alpha,double *A,int *lda,double *x,int *incx,
            double *beta,double *C,int *incy);
    void dgemm_(char* transa,char* transb,int *m,int *n,int *k,double *alpha,double *A,int *lda,
            double *B,int *ldb,double *beta,double *C,int *ldc);
    double ddot_(int *n,double *A,int *incx,double *B,int *incy);
    void dger_(int *m,int *n,double *alpha,double *A,int *incx,double *B,int *incy,double *C,int *lda);
    void dgetrf_(int *m,int *n,double *A,int *lda,int *ipivot,int *info);
    void dgetri_(int *n,double *A,int *lda,int *ipivot,double *work,int *lwork,int *info);
}

/** Prints a generic error message in cases when LAPACK reports an error. */
void lapack_fail() {
    fputs("LAPACK routine failed\n",stderr);
    exit(1);
}

/** Solves the matrix system Ax=b for the case when A is symmetric.
 * \param[in] n the dimension of the matrix.
 * \param[in] A the matrix terms; coefficients in the lower-triangular part
 *              will be ignored.
 * \param[in] x the source vector, which will be replaced by the solution
 *              when the routine exits.
 * \param[in] nrhs the number of right hand side vectors. */
void solve_sym_matrix(int n,double *A,double *x,int &nrhs) {

    // Create the temporary memory that LAPACK needs
    char uplo='u';
    int info,*ipiv=new int[n],lwork=64*n;
    double *work=new double[lwork];

    // Make LAPACK call
    dsysv_(&uplo,&n,&nrhs,A,&n,ipiv,x,&n,work,&lwork,&info);

    // Remove temporary memory
    delete [] work;
    delete [] ipiv;

    // Check for a non-zero value in info variable, indicating an error
    if(info!=0){
	    printf("info %d \n",info);
	    for(int j=0;j<n;j++) printf("Y[%d] %g\n",j,x[j]);
	    for(int j=0;j<n;j++)
		    for(int i=0;i<=j;i++)
			    printf("A[%d,%d] %g\n",i,j,A[j*n+i]);
	    lapack_fail();
    }
}

/* Factors a symmetric matrix.
 * \param[in] n the dimension of the matrix.
 * \param[in] A the matrix terms; coefficients in the lower-triangular part
 *              will be ignored. On exit the terms will contain the factorized
 *              matrix.
 * \param[in] ipiv the pivoting information of the factorization. */
void factor_sym(int n,double *A,int *ipiv) {

    // Create the temporary memory that LAPACK needs
    char uplo='u';
    int info,lwork=64*n;
    double *work=new double[lwork];

    // Perform LU decomposition
    dsytrf_(&uplo,&n,A,&n,ipiv,work,&lwork,&info);
    if(info!=0){
        printf("factor sym failed %d\n",info);
        lapack_fail();
    }

    // Remove temporary memory
    delete [] work;
}

/** Solves a symmetric linear system using a previously factored matrix.
 * \param[in] n the dimension of the matrix.
 * \param[in] A the factored matrix terms.
 * \param[in] ipiv the pivoting information.
 * \param[in] x the source vector, which will be replaced by the solution
 *              when the routine exits. */
void factor_sym_solve(int n,double *A,int *ipiv,double *x) {
    char uplo='u';
    int info,nrhs=1;
    dsytrs_(&uplo,&n,&nrhs,A,&n,ipiv,x,&n,&info);

    // Check for a non-zero value in info variable, indicating an error
    if(info!=0) lapack_fail();
}

/** Perform matrix multiplication.
 * \param[in] A,B matrices to multiply.
 * \param[in] C the return matrix.
 * \param[in] am the number of rows of the matrix A.
 * \param[in] an the number of columns of the matrix A.
 * \param[in] bn the number of columns of the matrix B. */
void matrix_mul(double *A,double *B,double *C,int am,int an,int bn,bool sym){
    int incx=1,incy=1;
    double alpha=1.,beta=0.;
    char uplo='u',side='l';

    if(sym){
        // symmetric matrix-vector multiply
        if(bn==1){
            dsymv_(&uplo,&am,&alpha,A,&am,B,&incx,&beta,C,&incy);

        }
        // symmetric matrix-general matrix multiply
        else{
            dsymm_(&side,&uplo,&am,&bn,&alpha,A,&am,B,&an,&beta,C,&am);

        }
    }
    else{
        // vector-vector multiply
        // A is 1 by an, B is an by 1
        if(am==1&&bn==1){
            *C=ddot_(&an,A,&incx,B,&incy);
        }
        // A is am by 1, B is 1 by bn
        if(an==1){
	    // Make sure C is all zero
	       for(int i=0;i<am*bn;i++) C[i]=0.;
            dger_(&am,&bn,&alpha,A,&incx,B,&incy,C,&am);
        }

        // matrix-vector multiply
        // AB=C, A is am by an matrix, B is an by 1 vector
        char tn='n';
        char tt='t';
        if(am!=1&&bn==1){
            dgemv_(&tn,&am,&an,&alpha,A,&am,B,&incx,&beta,C,&incy);
        }
        // AB=C, A is 1 by an vector, B is an by bn matrix
        if(am==1&&bn!=1){
            dgemv_(&tt,&bn,&an,&alpha,B,&bn,A,&incx,&beta,C,&incy);
        }

        // matrix-matrix multiply
        if(am!=1&&an!=1&&bn!=1){
            dgemm_(&tn,&tn,&am,&bn,&an,&alpha,A,&am,B,&an,&beta,C,&am);

        }
    }

}

/** Perform matrix inversion.
 * \param[in] A the matrix terms; coefficients in the lower-triangular part
 *              will be ignored. On exit the terms will contain the factorized
 *              matrix.
 * \param[in] n the dimension of the matrix. */
void matrix_inv(double *A,int n){

    // Create the temporary memory that LAPACK needs
    char uplo='u';
    int info,lwork=64*n;
    int *ipiv=new int[n];
    double *work=new double[lwork];

    // Perform LU decomposition
    dsytrf_(&uplo,&n,A,&n,ipiv,work,&lwork,&info);
    if(info!=0){
        puts("LU fail");
        lapack_fail();
    }

    // Perform inversion
    dsytri_(&uplo,&n,A,&n,ipiv,work,&info);
    if(info!=0){
        puts("Inversion fail");
        lapack_fail();
    }

    // Remove temporary memory
    delete [] ipiv;
    delete [] work;
}

/** Perform block matrix inversion
 * \param[in] A (n-1) by (n-1) inverse matrix.
 * \param[in] B n-1 vector.
 * \param[in] D scalar.
 * \param[in] n dimension of matrix.
 Return n by n inverse matrix. */
void matrix_update(double *A,double *B,double D,int n){
    double* E=new double[n-1];
    // Calculate E=A^-1 B
    matrix_mul(A,B,E,n-1,n-1,1,true);
    // for(int i=0;i<n-1;i++) printf("i %d E[i] %f B[i] %f\n",i,E[i],B[i]);
    // puts("a");
    double k;
    matrix_mul(E,B,&k,1,n-1,1,false);
    // puts("b");
    // printf("k %f\n",k);
    k=D-k;
    // printf("k %f\n",k);
    double* K1=new double[(n-1)*(n-1)];
    for(int i=0;i<(n-1)*(n-1);i++) K1[i]=0.;
    matrix_mul(E,E,K1,n-1,1,n-1,false);
    for(int j=0;j<n-1;j++){
        for(int i=0;i<=j;i++){
             // printf("i %d j %d K1[ij] %f\n",i,j,K1[j*(n-1)+i]);

        }
    }
     // puts("c");
    // Create a temporary A to store the new value
    double *AA=new double[n*n];
    // for(int i=0;i<n*n;i++) AA[i]=0.;
    for(int j=0;j<n-1;j++){
        for(int i=0;i<=j;i++){
            AA[j*n+i]=A[j*(n-1)+i]+K1[j*(n-1)+i]/k;
        }
    }
    for(int i=0;i<n-1;i++){
        AA[(n-1)*n+i]=-E[i]/k;
    }
    AA[(n-1)*n+n-1]=1./k;
    memcpy(A,AA,n*n*sizeof(double));

    delete [] AA;
    delete [] K1;
    delete [] E;

}
