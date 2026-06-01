#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "gpr.hh"
#include "turb_fluid.hh"
#include "common.hh"

inline double rnd(double b) {
    return (b/RAND_MAX)*static_cast<double>(rand());
}

int main() {
    gpr g(10);
    g.init(100,50,0.5,5);

    double se=0;
    for(int k=0;k<65536;k++) {
        double r=rnd(0.5),t=rnd(5),
               err=g.KF_eval(r,t)-g.KF_integrate(r,t);
        se+=err*err;
    }
    printf("%g\n",sqrt(se/65536));

    g.output_table("KF.bin");

    for(int i=0;i<100;i++){
	    double x=0.005*i;
	    double kk=g.KF_eval(x,0.8);
	    printf("%g %g\n",x,kk);
    }

}
