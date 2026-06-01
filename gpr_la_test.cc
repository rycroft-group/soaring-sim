#include <cstdio>

#include "k_func.hh"
#include "gpr.hh"

int main() {

    kernel_func KF(64,1,1,100,100,2,10);

    int m=32;

    gpr g(m,KF);

    for(int l=0;l<100000;l++) g.update_measurement(0,0.01*l,0,0,1,1,1);

    g.diag();
}
