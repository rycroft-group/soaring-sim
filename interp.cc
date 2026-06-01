#include <cstdio>

double f[5]={3.,1.,2.,4.,5.};

double state_revert(int level){
    switch(level){
        case 0: return 0.1;
        case 1: return 10.;
        case 2: return 50.;
        case 3: return 80.;
        default: return 100.;
    }
}

void state_level(double s,int &s0, int &s1){
    if(s<10){s0=0;s1=1;}
    else if(s<50){s0=1;s1=2;}
    else if(s<80){s0=2;s1=3;}
    else {s0=3;s1=4;}
}

/** Return value function V(r,s) */
double V(double s){
    int s0,s1;
    state_level(s,s0,s1);
    double f0=f[s0];
    double f1=f[s1];
    double w=(s-state_revert(s0))/(state_revert(s1)-state_revert(s0));
    return (1-w)*f0+w*f1;
}

int main() {
    for(double s=-10;s<120;s+=0.05)
        printf("%g %g\n",s,V(s));
}
