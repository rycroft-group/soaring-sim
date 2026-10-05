#include <cstdio>
#include <cstdlib>

#include "fileinfo.hh"
#include "soaring_sim.hh"
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc,char **argv) {

    // Check for a single argument for the filename of the configuration
    // file
    if(argc!=2) {
        fputs("./soar <config_file>\n",stderr);
        return 1;
    }

    // Set up multi-threaded FFTW computations, if available
#ifdef FFTW_OMP
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
#endif

    // Create the soaring_sim class, and read in all parameters from the
    // configuration file
    soaring_sim ss(argv[1]);

    // Choose a timestep, and print diagnostic information about the
    // simulation run
    ss.select_timestep(true);
    ss.print_info();

    // Run the simulation
    ss.run();
}
