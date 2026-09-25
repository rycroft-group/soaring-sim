#include "soaring_sim.hh"

#include <sys/types.h>
#include <sys/stat.h>
#include <limits>

/** Performs any global setup needed for file output. */
void soaring_sim::output_start() {

    // Create output directory
    mkdir(filename,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Set up memory for climb statistics
    if(fflags&1) {
        ginit=new double[8*gpt];
        cstats=new cli_stats[4*snaps+4];
        for(cli_stats *cs=cstats;cs<cstats+(4*snaps+4);cs++) cs->init();
    }

    // Set up memory for energy components
    if(fflags&2) {
        enstats=new cli_stats[3*snaps+3];
        for(cli_stats *en=enstats;en<enstats+(3*snaps+3);en++) en->init();
    }

    // Open file for GPR checksums
    if(fflags&4) {
        sprintf(fbuf,"%s/checksums",filename);
        checksums_file=safe_fopen(fbuf,"w");
    }
}

/** Performs any required setup for saving snapshots for a given fluid
 * simulation instance. */
void soaring_sim::snapshots_start() {

    // If needed, store the initial heights and energies of the gliders
    if(fflags&1) {
        glider_velocities();
        for(int i=0;i<gpt;i++) {
            ginit[2*i]=g[i].rz;
            ginit[2*i+1]=g[i].energy(wnd+3*i);
        }
    }

    // Print diagnostic message to checksums file if it's in use
    if(fflags&4) {
        fprintf(checksums_file,"# Wind field %d\n",fsim);
    }

    // If any snapshot types have been enabled that make multiple files, then
    // create a directory for storing them
    if(fflags&136) {
        sprintf(fbuf,"%s/f.%d",filename,fsim);
        mkdir(filename,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
    }

    // Open the glider position file if it's in use
    if(fflags&16) {
        sprintf(fbuf,"%s/glider_xyz.%d",filename,fsim);
        glider_xyz=safe_fopen(fbuf,"w");
    }

    // Open the minimal binary file if it's in use
    if(fflags&32) {
        sprintf(fbuf,"%s/glider_mb.%d",filename,fsim);
        glider_mb=safe_fopen(fbuf,"wb");
        fwrite(&gpt,sizeof(int),1,glider_mb);
        fwrite(&out_dur,sizeof(double),1,glider_mb);
    }

    // Open the larger binary file if it's in use
    if(fflags&64) {
        sprintf(fbuf,"%s/glider_lb.%d",filename,fsim);
        glider_lb=safe_fopen(fbuf,"wb");
        fwrite(&gpt,sizeof(int),1,glider_lb);
        fwrite(&out_dur,sizeof(double),1,glider_lb);
    }

    // Open the MCTS path file if it's in use
    if(fflags&1024) {
        sprintf(fbuf,"%s/path_info.%d",filename,fsim);
        path_file=safe_fopen(fbuf,"wb");
        double p_dur=c_dur*path_interval;
        fwrite(&gpt,sizeof(int),1,path_file);
        fwrite(&n_mcts,sizeof(int),1,path_file);
        fwrite(&mcts_depth,sizeof(int),1,path_file);
        fwrite(&out_pc,sizeof(int),1,path_file);
        fwrite(&p_dur,sizeof(double),1,path_file);
        fwrite(&out_dur,sizeof(double),1,path_file);
    }
}

/** Writes any selected snapshot files, and stores any information required for
 * the summary files.
 * \param[in] s the snapshot number. */
void soaring_sim::write_snapshots(int s) {
    double *dir=ginit+2*gpt,*dwnd=dir+3*gpt;

    // Store the information for outputting the glider climb summary files
    if(fflags&3) {

        // Calculate and store glider statistics: height, energy, and vertical
        // velocity
        for(int j=0;j<gpt;j++) {
            g[j].copy_pos(pos+3*j);
            g[j].copy_vel(dir+3*j);
        }
        tf->vel_dot_multi(gpt,pos,dir,wnd,dwnd);
        flag_current|=1;

        // Store the main glider climb summary file
        if(fflags&1) {
            cli_stats *cs=cstats+4*s;
            for(int i=0;i<gpt;i++) {

                // Use the cli_stats class to store the first and second moments,
                // and the extremal values
                cs->contrib(g[i].rz-ginit[2*i]);
                cs[1].contrib(g[i].energy(wnd+3*i)-ginit[2*i+1]);
                cs[2].contrib(g[i].uz);
                cs[3].contrib(gm->d_energy(g[i],wnd+3*i,dwnd+3*i));
            }
        }

        // Store the energy contributions
        if(fflags&2) {
            cli_stats *en=enstats+3*s;
            double alt,wacc;
            for(int i=0;i<gpt;i++) {

                // Use the cli_stats class to store the first and second moments,
                // and the extremal values
                en->contrib(gm->d_energy_component(g[i],wnd+3*i,dwnd+3*i,alt,wacc));
                en[1].contrib(alt);
                en[2].contrib(wacc);
            }
        }
    }

    // Output the wind cross-sections
    if(fflags&8) {
        for(unsigned int o=0;o<w_cs.size();o++) {
            int ind=w_cs[o]>>4,dir=w_cs[o]&3,cmp=(w_cs[o]>>2)&3;
            sprintf(fbuf,"%s/f.%d/w%c_%c.%d.%d",
                    filename,fsim,'x'+cmp,'x'+dir,ind,s);
            switch(dir) {
                case 0: tf->output_x(fbuf,cmp,ind);break;
                case 1: tf->output_y(fbuf,cmp,ind);break;
                case 2: tf->output_z(fbuf,cmp,ind);break;
                default: break;
            }
        }
    }

    // Output the glider positions as a text file
    if(fflags&16) {
        fprintf(glider_xyz,"%d %g",s,out_dur*s);
        for(int i=0;i<gpt;i++) g[i].write_pos(glider_xyz);
        fputc('\n',glider_xyz);
        fflush(glider_xyz);
    }

    // Output the glider positions as a minimal binary file
    if(fflags&32)
        for(int i=0;i<gpt;i++) g[i].write_min_binary(glider_mb);

    // Output the glider positions and wind data as a larger binary file
    if(fflags&64) {
        glider_velocities();
        for(int i=0;i<gpt;i++) g[i].write_binary(glider_lb,wnd+3*i);
    }

    // Output the complete wind state
    if(fflags&128) {
        sprintf(fbuf,"%s/f.%d/wmodes.%d",filename,fsim,s);
        tf->save(fbuf);
    }
}

/** Performs any required finalization for completing snapshots for a given
 * fluid instance. */
void soaring_sim::snapshots_end() {
    if(fflags&4) {
        for(int j=0;j<gpt;j++) mem[j]->diagnostics(checksums_file);
        fputc('\n',checksums_file);
    }
    if(fflags&16) fclose(glider_xyz);
    if(fflags&32) fclose(glider_mb);
    if(fflags&64) fclose(glider_lb);
    if(fflags&1024) fclose(path_file);
}

/** Writes any summary files, and performs any finalization for output at the
 * end of the simulation. */
void soaring_sim::output_end() {

    // Output the glider climb summary
    if(fflags&1) {
        sprintf(fbuf,"%s/cstats",filename);
        FILE *fp=safe_fopen(fbuf,"w");

        // Use the first and second moments to compute the mean and standard
        // deviation the glider climb as a function of time
        cli_stats *cs=cstats;
        for(int i=0;i<=snaps;i++,cs+=4) {
            double nor=1./(num_trials*gpt),mu0,sig0,mu1,sig1,mu2,sig2,mu3,sig3;
            cs->mu_sig(nor,mu0,sig0);
            cs[1].mu_sig(nor,mu1,sig1);
            cs[2].mu_sig(nor,mu2,sig2);
            cs[3].mu_sig(nor,mu3,sig3);
            fprintf(fp,"%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                    i,out_dur*i,mu0,sig0,cs->mi,cs->ma,
                    mu1,sig1,cs[1].mi,cs[1].ma,
                    mu2,sig2,cs[2].mi,cs[2].ma,
                    mu3,sig3,cs[3].mi,cs[3].ma);
        }
        fclose(fp);

        // Deallocate temporary memory
        delete [] cstats;
        delete [] ginit;
    }

    // Output the energy components
    if(fflags&2) {
        sprintf(fbuf,"%s/en_components",filename);
        FILE *fp=safe_fopen(fbuf,"w");
        cli_stats *en=enstats;
        for(int i=0;i<=snaps;i++,en+=3) {
            double nor=1./(num_trials*gpt),mu0,sig0,mu1,sig1,mu2,sig2;
            en->mu_sig(nor,mu0,sig0);
            en[1].mu_sig(nor,mu1,sig1);
            en[2].mu_sig(nor,mu2,sig2);
            fprintf(fp,"%d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                    i,out_dur*i,mu0,sig0,en->mi,en->ma,
                    mu1,sig1,en[1].mi,en[1].ma,
                    mu2,sig2,en[2].mi,en[2].ma);
        }
        fclose(fp);

        // Deallocate temporary memory
        delete [] enstats;
    }

    // Close the GPR checksums file
    if(fflags&4) {
        fprintf(checksums_file,"# Kernel table misses: %d\n",KF->miss);
        fclose(checksums_file);
    }

    // Output the wind correlation information
    if(fflags&256) {
        sprintf(fbuf,"%s/correl",filename);
        FILE *fp=safe_fopen(fbuf,"w");
        wic->output(fp);
        fclose(fp);
    }

    // Output the MCTS tree info
    if(fflags&512) {
        sprintf(fbuf,"%s/tree_info",filename);
        FILE *fp=safe_fopen(fbuf,"w");

        // Loop over the range of MCTS depths
        double mu,sig;
        mti_stats ms;
        for(int o=0;o<=mcts_depth;o++) {

            // At each level in the tree, aggregate all of the accumulators
            // across the MCTS classes
            ms.init();
            for(int i=0;i<mctst;i++) ms.add(mc[i]->ms[o]);
            ms.mu_sig(mu,sig);

            // Print the statistics
            fprintf(fp,"%d %ld %g %g %d %d\n",o,ms.n,mu,sig,ms.mi,ms.ma);
        }
        fclose(fp);
    }
}
