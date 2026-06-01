#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "gpr.hh"
#include "common.hh"

/** Calculate uncertainty of GPR for helices with differencet curvature. */

inline double rnd(double b) {
    return (b/RAND_MAX)*static_cast<double>(rand());
}

int main() {
    gpr g(100);
    g.init(200,50,sqrt(3.),5);

    // Create samples
    int num_x=50;
    int num_xs=1;
    measure_info* X=new measure_info[num_x];
    measure_info* X_star=new measure_info[num_xs];

    const double a=0.05;
    const double b=0.01;
    const double ds=0.1;
    // Want ds=dt*sqrt(a*a+b*b) constant
    double dt=ds/sqrt(a*a+b*b);
    double curvature=a/(a*a+b*b);
    printf("curvature %g torsion %g dt %g duration %g\n",curvature,b/(a*a+b*b),dt,dt*num_x/(2*M_PI));

    double x,y,z,t=0.;
    for(int i=0;i<num_x;i++){
        //x=a*cos(t);
	x=0.01*i;
        // y=a*sin(t);
       y=0.0*i;
        //z=b*t;
       z=0.;
        // Want dt*sqrt(a*a+b*b) constant
        printf("%4g %4g %4g\n",x,y,z);
        X[i].set(x,y,z,0.,0.,0.,0.);
        t+=dt;

    }

    for(int i=0;i<num_xs;i++){
        //x=a*cos(t);
        x=1.;
	// y=a*sin(t);
        y=0.;
         //z=b*t;
        z=0.;
        printf("%4g %4g %4g\n",x,y,z);
        X_star[i].set(x,y,z,0.,0.,0.,0.);
        t+=dt;
    }
    double* result=new double[num_xs*num_xs];
    g.predict_var(num_x,num_xs,X,X_star,result);

    for(int i=0;i<num_xs*num_xs;i++) printf("result %g\n",result[i]);

    delete [] X;
    delete [] X_star;

}
