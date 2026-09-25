#ifndef COMMON_HH
#define COMMON_HH

#include <cstdio>

// Set up a timing routine, using either OpenMP or the standard C++ library
#ifdef _OPENMP
#include "omp.h"
inline double wtime() {return omp_get_wtime();}
#else
#include <ctime>
inline double wtime() {return double(clock())/CLOCKS_PER_SEC;}
#endif

/** The numerical integration type. */
enum integration_type {
    it_unset, it_euler, it_improv_e
};

FILE* safe_fopen(const char* filename,const char* mode);
void fatal_error(const char *p,int code);

#endif
