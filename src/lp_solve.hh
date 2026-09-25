#ifndef LP_SOLVE_HH
#define LP_SOLVE_HH

void lapack_fail();
void solve_sym_matrix(int n,double *A,double *x,int &nrhs);
void factor_sym(int n,double *A,int *ipiv);
void factor_sym_solve(int n,double *A,int *ipiv,double *x);

void matrix_mul(double *A,double *B,double *C,int am,int an,int bn,bool sym);
void matrix_inv(double *A,int n);
void matrix_update(double *A,double *B,double D,int n);
#endif
