#include "game.hh"
#include "game_mcts.hh"
#include "turb_fluid.hh"
#include "gpr.hh"
#include <gsl/gsl_rng.h>
#include "common.hh"
#include <cstring>
#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _OPENMP
#include "omp.h"
#endif

int main(){
    gsl_rng *rng;
    rng=gsl_rng_alloc(gsl_rng_taus2);

    int action_size=3;
#pragma omp parallel for
    for(int l=0;l<20;l++){
        gpr g(20);
        g.init();
        turb_fluid_grid tf(256,256,256,0,50,0,50,0,50,1.,10.);
        tf.init_steady_state();

        // Initialize glider
        game_glider g0(1.,1/15.,0.05,0,true,tf,g);
        game_glider g1(1.,1/15.,0.05,0,true,tf,g);
        game_glider g2(1.,1/15.,0.05,0,true,tf,g);
        game_glider g3(1.,1/15.,0.05,0,true,tf,g);
        game_glider g_random(1.,1/15.,0.05,0,true,tf,g);
        game_glider g_policy(1.,1/15.,0.05,0,true,tf,g);

        // Start at random locations
        double rx_,ry_,rz_;//,ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_,bank_;
        rx_=50.*gsl_rng_uniform(rng);
        ry_=50.*gsl_rng_uniform(rng);
        rz_=50.*gsl_rng_uniform(rng);
        //rz_=0.;
        g0.init(rx_,ry_,rz_,1.,0.,-1./15);
        g1.init(rx_,ry_,rz_,1.,0.,-1./15);
        g2.init(rx_,ry_,rz_,1.,0.,-1./15);
        g3.init(rx_,ry_,rz_,1.,0.,-1./15);
        g_random.init(rx_,ry_,rz_,1.,0.,-1./15);
        g_policy.init(rx_,ry_,rz_,1.,0.,-1./15);

        char fn[50];
        char buffer[100];
        // Open file and store initial positions
        sprintf(fn,"test_random/tf_result_%d",l);
        mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
        sprintf(buffer,"%s/g0.dat",fn);
        FILE *fp0=safe_fopen(buffer,"w");
        sprintf(buffer,"%s/g1.dat",fn);
        FILE *fp1=safe_fopen(buffer,"w");
        sprintf(buffer,"%s/g2.dat",fn);
        FILE *fp2=safe_fopen(buffer,"w");
        sprintf(buffer,"%s/g3.dat",fn);
        FILE *fp3=safe_fopen(buffer,"w");
        sprintf(buffer,"%s/g_random.dat",fn);
        FILE *fp4=safe_fopen(buffer,"w");
        sprintf(buffer,"%s/g_policy.dat",fn);
        FILE *fp5=safe_fopen(buffer,"w");

        int i=1,action;
        int action_scale=50;
        g0.bank=0;
        g1.bank=g1.mu_f(2);
        g2.bank=g2.mu_f(2)*2.;
        g3.bank=g3.mu_f(2)*3.;
        while(i<5000){
            if(i%1==0) printf("time frame %d\n",i);
            action=g_policy.pick_action_policy(i);
            g_policy.bank+=g_policy.mu_f(action);
            for(int k=0;k<action_scale;k++,i++){
                 g_policy.play(action,1,i);
                 g_policy.output(i+1,action,fp5);
            }
            action=gsl_rng_uniform_int(rng,action_size);
            g_random.bank+=g_random.mu_f(action);
            i-=action_scale;
            for(int k=0;k<action_scale;k++,i++){
                g0.play(action,1,i);
                g0.output(i+1,action,fp0);
                g1.play(action,1,i);
                            g1.output(i+1,action,fp1);
                g2.play(action,1,i);
                            g2.output(i+1,action,fp2);
                g3.play(action,1,i);
                            g3.output(i+1,action,fp3);
                g_random.play(action,1,i);
                            g_random.output(i+1,action,fp4);
            }

        }
        if(fp0!=NULL) fclose(fp0);
        if(fp1!=NULL) fclose(fp1);
        if(fp2!=NULL) fclose(fp2);
        if(fp3!=NULL) fclose(fp3);
        if(fp4!=NULL) fclose(fp4);
        if(fp5!=NULL) fclose(fp5);
    }
}
