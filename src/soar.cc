#include <cstdio>
#include <cstdlib>

#include "fileinfo.hh"
#include "soaring_sim.hh"

int main(int argc,char **argv) {

    if(argc!=2) {
        fputs("./fi_test <config_file>\n",stderr);
        return 1;
    }
    soaring_sim ss(argv[1]);

    ss.select_timestep(true);
    ss.print_info();

    //ss.tf->init_steady_state();
    //ss.tf->transform();

    //double ubar,vbar,wbar,urms,vrms,wrms;
    //ss.tf->grid_vel_stats(ubar,vbar,wbar,urms,vrms,wrms);
    //printf("%g %g %g %g\n",urms,vrms,wrms,sqrt(urms*urms+vrms*vrms+wrms*wrms));
    ss.run();
}
