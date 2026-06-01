#include "gpr.hh"
#include "lp_solve.hh"
#include "common.hh"

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// Tell the compiler about the existence of the required LAPACK functions
extern "C" {
    void dsytrf_(char* uplo,int *n,double *a,int *lda,int *ipivot,double *work,
            int *lwork,int *info);
    void dsytri_(char* uplo,int *n,double *A,int *lda,int *ipivot,double *work,int *info);
    void dsymm_(char* side,char* uplo,int *m,int *n,double *alpha,double* A,int *lda,
            double *B,int *ldb,double *beta,double *C,int *ldc);
}

/** The Gaussian process regression constructor sets up class constants and
 * dynamically allocates memory.
 * \param[in] m_ The maximum number of measurements to use in the regression.
 * \param[in] KF_ a reference to the kernel function evaluation class.
 * \param[in] ker_tol the tolerance to use on detecting ill-conditioned covariance
 *                     updates to skip.
 * \param[in] ker_step a typical distance covered between a pair of measurements,
 *                     used to estimate the kernel decay.
 * \param[in] full_compute_ Whether to always bypass the Woodbury formula and do
 *                          a full LAPACK computation each time during updates. */
gpr::gpr(int m_,kernel_func &KF_,double ker_tol,double ker_step,bool full_compute_)
    : m(m_), mm(m*m), v(0), lwork(64*m), full_compute(full_compute_),
    full(false), M(new measure_info[m]), A(new double[mm]), B(new double[mm]),
    W(new double[3*m]), X(new double[3*m]), ipiv(new int[m]),
    work(new double[lwork]), KF(KF_), ker_det_tol(ker_tol),
    ker_det_tol2(ker_tol*(2-2*KF.eval(ker_step,0))) {

    // The diagonal entries of the A matrix are always 1 due to
    // self-correlation. Set them here once.
    for(int i=0;i<mm;i+=m+1) A[i]=1;

    // Clear counters
    *co=co[1]=co[2]=co[3]=0;
}

/** The class destructor frees the dynamically allocated memory. */
gpr::~gpr() {
    delete [] work;
    delete [] ipiv;
    delete [] X;
    delete [] W;
    delete [] B;
    delete [] A;
    delete [] M;
}

/** Adds a wind measurement to the Gaussian process regression.
 * \param[in] (x,y,z) the position of the measurement.
 * \param[in] t the time of the measurement.
 * \param[in] (wx,wy,wz) the measured wind vector.*/
void gpr::add_measurement(double x,double y,double z,double t,double wx,double wy,double wz) {

    /// Store the measurement position/time and the wind information
    M[v].set(x,y,z,t);
    W[v]=wx;W[m+v]=wy;W[2*m+v]=wz;

    // Update entries of the covariance matrix
    for(int i=0;i<v;i++) A[v*m+i]=K(M[i],x,y,z,t);
    if(full) for(int i=v+1;i<m;i++) A[i*m+v]=K(M[i],x,y,z,t);

    // Increment the row, looping back to zero if the matrix is full
    if(++v==m) {full=true;v=0;}
}

/** Computes the inverse of the currently used portion of the covariance matrix.
 * \param[in] Q the place to store the inverse. */
void gpr::compute_inverse(double *Q) {

    // Copy A matrix into A inverse
    int info,n=full?m:v;
    for(int j=0;j<n;j++) {
        for(int i=0;i<j;i++) Q[j*m+i]=A[j*m+i];
        Q[j*(m+1)]=1;
    }

    // Perform LU decomposition
    char uplo='u';
    dsytrf_(&uplo,&n,Q,&m,ipiv,work,&lwork,&info);
    if(info!=0) fatal_error("LU calculation failed",1);

    // Perform inversion
    dsytri_(&uplo,&n,Q,&m,ipiv,work,&info);
    if(info!=0) fatal_error("Matrix inversion failed",1);
}

/** Calculates the X array that is used for making wind predictions. */
void gpr::compute_X() {
    double alpha=1,beta=0;
    char uplo='u',side='l';
    int nrhs=3,n=full?m:v;
    dsymm_(&side,&uplo,&n,&nrhs,&alpha,B,&m,W,&m,&beta,X,&m);
}

/** Adds a wind measurement to the Gaussian process regression, and updates the
 * kernel calculations ready for performing predictions.
 * \param[in] (x,y,z) the position of the measurement.
 * \param[in] t the time of the measurement.
 * \param[in] (wx,wy,wz) the measured wind vector.*/
void gpr::update_measurement(double x,double y,double z,double t,double wx,double wy,double wz) {
    int n;

    /// Store the measurement position/time and the wind information
    M[v].set(x,y,z,t);
    W[v]=wx;W[m+v]=wy;W[2*m+v]=wz;

    if(full) {

        // Set up pointers to temporary workspace
        double *R=work,*S=R+m,*V=S+m,KK,alpha,beta,gamma=0,sigma;

        // Extract row from the A^{-1} matrix
        for(int i=0;i<v;i++) {
            KK=K(M[i],x,y,z,t);
            V[i]=KK-A[i+v*m];
            A[i+v*m]=KK;
        }
        alpha=B[v+v*m];V[v]=0;
        for(int i=v+1;i<m;i++) {
            KK=K(M[i],x,y,z,t);
            V[i]=KK-A[v+i*m];
            A[v+i*m]=KK;
        }

        // Compute the A^{-1} V matrix vector product
        for(int i=0;i<m;i++) {
            S[i]=0.;
            for(int j=0;j<m;j++) S[i]+=(j<i?B[j+i*m]:B[i+j*m])*V[j];
            gamma+=V[i]*S[i];
        }

        // Compute the determinant of the 2x2 matrix in the Woodbury formula.
        // If this is small or negative, then the A matrix becomes
        // ill-conditioned or loses positive definiteness.
        beta=1+S[v];sigma=beta*beta-alpha*gamma;
        if(sigma<ker_det_tol) {
            co[0]++;
            M[v].nullify();

            // Zero out the off-diagonal entries of the row+column
            for(int i=0;i<v;i++) A[i+v*m]=0;
            for(int i=v+1;i<m;i++) A[v+i*m]=0;
            compute_inverse(B);

        } else {
            co[1]++;
            if(full_compute) compute_inverse(B);
            else {
                for(int i=0;i<v;i++) R[i]=B[i+v*m];
                R[v]=alpha;
                for(int i=v+1;i<m;i++) R[i]=B[v+i*m];

                // Apply the Woodbury formula
                sigma=1./sigma;
                gamma*=sigma;beta*=sigma;alpha*=sigma;
                for(int j=0;j<m;j++) for(int i=0;i<=j;i++)
                    B[j*m+i]+=gamma*R[i]*R[j]-beta*(R[i]*S[j]+S[i]*R[j])+alpha*S[i]*S[j];
            }
        }

        if(++v==m) v=0;
        n=m;
    } else {
        for(int i=0;i<v;i++) A[v*m+i]=K(M[i],x,y,z,t);

        // Compute the A^{-1} V matrix vector product
        double gamma=0,*S=work;
        for(int i=0;i<v;i++) {
            S[i]=0.;
            for(int j=0;j<v;j++) S[i]+=(j<i?B[j+i*m]:B[i+j*m])*A[v*m+j];
            gamma+=S[i]*A[v*m+i];
        }

        double sigma=1-gamma;
        if(sigma<ker_det_tol2) {

            // Add an identity line to the matrices, corresponding to a skipped
            // measurement
            for(int i=0;i<v;i++) A[v*m+i]=B[v*m+i]=0;
            B[v*m+v]=1;
            co[2]++;
            n=++v;
        } else {

            // Apply the Woodbury formula for this restricted case
            if(full_compute) {
                n=++v;
                compute_inverse(B);
            } else {
                sigma=1./sigma;
                for(int j=0;j<v;j++) {
                    for(int i=0;i<=j;i++) B[j*m+i]+=sigma*S[i]*S[j];
                    B[v*m+j]=-sigma*S[j];
                }
                B[v*m+v]=1+gamma*sigma;
                n=++v;
            }
            co[3]++;
        }
        if(n==m) {full=true;v=0;}
    }

    // Perform the matrix inverse multiplication with the RHS
    double alpha=1,beta=0;
    int nrhs=3;
    char uplo='u',side='l';
    dsymm_(&side,&uplo,&n,&nrhs,&alpha,B,&m,W,&m,&beta,X,&m);
}

/** Predicts the wind from the given measurements.
 * \param[in] (x,y,z) the position of the measurement.
 * \param[in] t the time of the measurement.
 * \param[out] (wx,wy,wz) the predicted wind vector. */
void gpr::predict(double x,double y,double z,double t,double &wx,double &wy,double &wz) {

    // Assemble the predicted wind velocity
    wx=wy=wz=0;
    int n=full?m:v;
    for(int i=0;i<n;i++) {
        double k=K(M[i],x,y,z,t);
        wx+=X[i]*k;
        wy+=X[i+m]*k;
        wz+=X[i+2*m]*k;
    }
}

void gpr::print(double *Q) {
    int n=full?m:v;
    for(int i=0;i<n;i++) {
        printf("[");
        for(int j=0;j<n;j++) {
            if(i>j) printf(" *******");
            else printf(" %7.3g",Q[j*m+i]);
        }
        puts(" ]");
    }
}

double gpr::checksum() {
    double *C=new double[mm];
    compute_inverse(C);

    double s=0;
    for(int j=0;j<m;j++) for(int i=0;i<=j;i++) {
        C[i+m*j]-=B[i+m*j];
        s+=C[i+m*j]*C[i+m*j];
    }

    delete [] C;
    return 2.*s/(m*(m+1));
}

void gpr::diagnostics(FILE *fp) {
    fprintf(fp,"%.14g %ld %ld %ld %ld\n",checksum(),*co,co[1],co[2],co[3]);
    *co=co[1]=co[2]=co[3]=0;
}
