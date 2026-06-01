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

int main(int argc, char **argv)
{
	if(argc!=3){
        fputs("Syntax: ./game_frozen <integer> <float> \n",stderr);
        return 1;
    }

    int memory=atoi(argv[1]);
    double time_scale=atof(argv[2]);

    // Perform some basic checks on the numbers
    if(memory<0) {
        fputs("<integer> must be positive\n",stderr);
        return 1;
    }
    if(time_scale<0) {
        fputs("<float> must be positive\n",stderr);
        return 1;
    }

	double t0=wtime();
	gsl_rng *rng;
	rng=gsl_rng_alloc(gsl_rng_taus2);

	int action_size=3;
	int wind_state=0; // 0: partial information; 1: stable currents; 2: full information.
	int depth=1000;
	int nmode=64;
	const int ngl=10;
	int action_scale=50;
	const double box_size=50.;
	double alpha=M_PI*M_PI/(9.*(pow(sqrt(3.)*M_PI,-2/3.)-pow(sqrt(3.)*M_PI*nmode,-2/3.)));

	bool frozen=true;
	double wind_val=0.05;
	double *r=new double[3*ngl];
	double *wind=new double[3*ngl];
	char fn[100];
	char buffer[150];
	bool debug=false;

	for(int l=0;l<20;l++){
		printf("tf %d m %d d %d frozen %d wind_val %g time scale %g \n",l,memory,depth,frozen,wind_val,time_scale);
		turb_fluid_grid tf(nmode,nmode,nmode,0,box_size,0,box_size,0,box_size,0.67,alpha);
		tf.init_steady_state();
		tf.transform();

		gpr **gp_test=new gpr*[ngl];
		game_glider **gl_frozen=new game_glider*[ngl];
		FILE **fp1=new FILE*[ngl];
		FILE **fp2=new FILE*[ngl];

		// Initialization
		for(int o=0;o<ngl;o++){
			gp_test[o]=new gpr(memory,time_scale);
			gp_test[o]->init();
			gl_frozen[o]=new game_glider(1.,1/15.,wind_val,wind_state,frozen,tf,*gp_test[o]);
			double rx_=box_size*gsl_rng_uniform(rng);
  			double ry_=box_size*gsl_rng_uniform(rng);
            double rz_=box_size*gsl_rng_uniform(rng);
            gl_frozen[o]->init(rx_,ry_,rz_,1.,0.,-1./15);

            // Open files
            sprintf(fn,"test/frozen/m%d_d%d_a%d/tf_result_%d",memory,depth,action_scale,l*ngl+o);
            mkdir(fn,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
            sprintf(buffer,"%s/game_frozen.dat",fn);
            fp1[o]=safe_fopen(buffer,"w");
	        sprintf(buffer,"%s/game_frozen_sim.dat",fn);
            fp2[o]=safe_fopen(buffer,"w");
            gl_frozen[o]->output(0,0,fp1[o]);
	     // Get current wind status and add measurement
	}
	    for(int o=0;o<ngl;o++){
                r[3*o]=gl_frozen[o]->rx;
                r[3*o+1]=gl_frozen[o]->ry;
                r[3*o+2]=gl_frozen[o]->rz;
            }
            // Calculate the velocities at the random sample points
            tf.allocate_vel_table(ngl);
            tf.vel_multi(ngl,r,wind);
            for(int o=0;o<ngl;o++){
                gl_frozen[o]->wx=wind_val*wind[3*o];
                gl_frozen[o]->wy=wind_val*wind[3*o+1];
                gl_frozen[o]->wz=wind_val*wind[3*o+2];
                gp_test[o]->add_measurement(r[3*o]/box_size,r[3*o+1]/box_size,r[3*o+2]/box_size,0.,wind_val*wind[3*o],wind_val*wind[3*o+1],wind_val*wind[3*o+2]);

            }

		int i=1;
		// Start simulation
		while(i<50000){
			tf.transform();
			#pragma omp parallel for schedule(dynamic)
			// MCTS simulation
			for(int o=0;o<ngl;o++){
				double rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,wx_,wy_,wz_,bank_;
				game_mcts mcts(action_size,depth,*gl_frozen[o]);
				mcts.reset();

				// Save the current status
				rx_=gl_frozen[o]->rx,ry_=gl_frozen[o]->ry,rz_=gl_frozen[o]->rz;
                ux_=gl_frozen[o]->ux,uy_=gl_frozen[o]->uy,uz_=gl_frozen[o]->uz;
                vx_=gl_frozen[o]->vx,vy_=gl_frozen[o]->vy,vz_=gl_frozen[o]->vz;
                wx_=gl_frozen[o]->wx,wy_=gl_frozen[o]->wy,wz_=gl_frozen[o]->wz;
                bank_=gl_frozen[o]->bank;

                // Start searching
                bool out=false;
                for(int j=0;j<10000;j++){
			out=false;
			if((j%1000==0)&&(i%1000==1)) out=true;
                	mcts.simulate_path_glider(1,i,out,fp2[o]);
                }
                // Find the best action
                int action=mcts.pick_best(0);
                //Restore the current status
                gl_frozen[o]->reset(rx_,ry_,rz_,ux_,uy_,uz_,vx_,vy_,vz_,bank_);
                // Take an action
                gl_frozen[o]->bank+=gl_frozen[o]->mu_f(action);

			}

			// Play
			for(int k=0;k<action_scale;k++,i++){
				// Get positions
				for(int o=0;o<ngl;o++){
					r[3*o]=gl_frozen[o]->rx;
                    r[3*o+1]=gl_frozen[o]->ry;
                    r[3*o+2]=gl_frozen[o]->rz;
				}
				// Calculate the velocities at the random sample points
                tf.allocate_vel_table(ngl);
                tf.vel_multi(ngl,r,wind);
                for(int o=0;o<ngl;o++){
                	// Set wind velocities
                	double wx_,wy_,wz_;
					// tf.vel(r[3*o],r[3*o+1],r[3*o+2],wx_,wy_,wz_);
					gl_frozen[o]->wx=wind_val*wind[3*o];
                    gl_frozen[o]->wy=wind_val*wind[3*o+1];
                    gl_frozen[o]->wz=wind_val*wind[3*o+2];
                    if(k==action_scale-1){
			    double time=0.;
			    if(!frozen) time=i*gl_frozen[o]->dt;
			    gp_test[o]->add_measurement(r[3*o]/box_size,r[3*o+1]/box_size,r[3*o+2]/box_size,time,wind_val*wind[3*o],wind_val*wind[3*o+1],wind_val*wind[3*o+2]);
		    }
		    gl_frozen[o]->wind=3;
                    gl_frozen[o]->play(0,1,i);
                    gl_frozen[o]->wind=wind_state;
                    // Store the move
                    gl_frozen[o]->output(i+1,0,fp1[o]);
                }
			}

		}

		// Delete pointers
		for(int o=0;o<ngl;o++){
			delete gl_frozen[o];
			delete gp_test[o];
			if(fp1[o]!=NULL) fclose(fp1[o]);
			if(fp2[o]!=NULL) fclose(fp2[o]);
		}

		delete [] gl_frozen;
		delete [] gp_test;
		delete [] fp1;
		delete [] fp2;

	}

	delete [] wind;
	delete [] r;
	double t1=wtime();
	printf("running time is %g\n",t1-t0);

}
