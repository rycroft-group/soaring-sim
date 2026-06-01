#include <cstdio>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.hh"
#include "k_func.hh"

// The duration to integrate over
const double duration=0.1;

// The number of instances
const int ins=81;

int main(int argc,char **argv) {

    //const int nbin=601;
    const double rmax=90;//h=rmax/(nbin-1);

    kernel_rt KF(64,50,1,100,100,rmax,duration*(ins-1));
    kernel_r KFr(64,50,100,rmax);

    KF.output_table("KF");
    KFr.output_table("KFr");
}
