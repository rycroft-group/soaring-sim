#include "fileinfo.hh"

#include <cmath>
#include <limits>

/** The class constructor reads parameters from an input file, converts some
 * physical units into simulation units, and performs a number of consistency
 * checks.
 * \param[in] infile the name of the input file to read. */
fileinfo::fileinfo(const char* infile) : num_trials(-1), gpt(-1), nx(-1),
    ny(-1), nz(-1), gprm(-1), gs_pc(-1), n_mcts(-1), mcts_depth(-1),
    mcts_ex_fac(-1), h_segs(-1), gpr_rs(-1), gpr_ts(-1), bypass_mr(false),
    gpr_full_compute(false), fflags(0), base_seed(1), l_phys(-1), t_phys(-1),
    v_phys(-1), g_phys(-1), nu_phys(-1), lx(-1), ly(-1), lz(-1), w_Cinv(-1),
    w_alpha(-1), w_rms(-1), gpr_dist_pad(-1), gpr_ker_tol(0), c_dur(-1),
    duration(-1), wmodel(wm_unset), itype(it_unset), ptype(pt_unset) {
    int br_min=1,br_max=0,bc_min=1,bc_max=0;
    double bstep=-1,c_D=-1,c_L=-1,c_dur_phys=-1,w_rms_phys=-1,
           wic_tcut_phys=-1,mti_tcut_phys=-1,w_Cinv_phys=-1,duration_phys=-1;

    // Check that the filename ends in '.cfg'
    int l=strlen(infile),ln=1;
    if(l<4) fatal_error("Filename is too short",1);
    const char* ip=infile+l-4;
    if(*ip!='.'||ip[1]!='c'||ip[2]!='f'||ip[3]!='g')
        fatal_error("Filename must end in '.cfg'",1);

    // Assemble output filename by replacing '.cfg' with '.odr'. Allocate
    // buffer for creating output filenames.
    filename=new char[l+1];
    fbuf=new char[l+fileinfo_fbuf_pad_size];
    memcpy(filename,infile,l-3);
    char *fp=filename+l-3;
    *fp='o';fp[1]='d';fp[2]='r';fp[3]=0;

    // Open the input file and read
    FILE *f=safe_fopen(infile,"r");
    char *buf=new char[fileinfo_buf_size],*bp;
    while(!feof(f)) {
        if(fgets(buf,fileinfo_buf_size,f)==NULL) break;

        // Locate comments and remove by replacing comment character
        // with a null character
        bp=buf;
        while((*bp)!=0) {
            if(*bp=='#') {*bp=0;break;}
            bp++;
        }

        // Separate entries in the file by tabs and spaces; if no entries are
        // available then skip this line
        bp=strtok(buf," \t\n");
        if(bp==NULL) {ln++;continue;}

        // Look for a known keyword, and read in any extra values
        if(se(bp,"num_trials")) {
            num_trials=final_int(ln);
            if(num_trials<=0) fatal_error("Number of trials must be positive",1);
        } else if(se(bp,"gliders_per_trial")) {
            gpt=final_int(ln);
            if(gpt<=0) fatal_error("Gliders per trial must be positive",1);
        } else if(se(bp,"wf_modes")) {
            nx=ny=nz=final_int(ln);
            if(nx<=3) fatal_error("Number of modes too small",1);
        } else if(se(bp,"wf_modes_xyz")) {
            nx=atoi(next_token(ln));
            ny=atoi(next_token(ln));
            nz=final_int(ln);
            if(nx<=3||ny<=3||nz<=3) fatal_error("Number of modes too small",1);
        } else if(se(bp,"bank_angle_range")) {
            br_min=atoi(next_token(ln));
            br_max=final_int(ln);
            if(br_min>br_max) fatal_error("Min bank angle larger than max",1);
        } else if(se(bp,"bank_angle_control")) {
            bc_min=atoi(next_token(ln));
            bc_max=final_int(ln);
            if(bc_min>bc_max) fatal_error("Min bank angle larger than max",1);
        } else if(se(bp,"gpr_memory")) {
            gprm=final_int(ln);
            if(gprm<=0) fatal_error("GPR memory must be positive",1);
        } else if(se(bp,"gpr_samp_per_ctl")) {
            gs_pc=final_int(ln);
            if(gs_pc<=0) fatal_error("GPR samples per control period must be positive",1);
        } else if(se(bp,"output_per_ctl")) {
            out_pc=final_int(ln);
            if(out_pc<=0) fatal_error("Outputs per control period must be positive",1);
        } else if(se(bp,"n_mcts")) {
            n_mcts=final_int(ln);
            if(n_mcts<=0) fatal_error("Number of MCTS trials must be positive",1);
        } else if(se(bp,"bypass_mr")) {
            bypass_mr=true;
            check_no_more(ln);
        } else if(se(bp,"gpr_full_compute")) {
            gpr_full_compute=true;
            check_no_more(ln);
        } else if(se(bp,"mcts_depth")) {
            mcts_depth=final_int(ln);
            if(mcts_depth<=0) fatal_error("MCTS depth must be positive",1);
        } else if(se(bp,"mcts_ex_fac")) {
            mcts_ex_fac=atof(next_token(ln));
            check_no_more(ln);
            if(mcts_ex_fac<=0) fatal_error("MCTS exploration factor must be positive",1);
        } else if(se(bp,"hermite_segments")) {
            h_segs=final_int(ln);
            if(h_segs<=0) fatal_error("The number of Hermite segments must be positive",1);
        } else if(se(bp,"random_seed")) {
            base_seed=atol(next_token(ln));
            check_no_more(ln);
        } else if(se(bp,"l_phys")) l_phys=final_double(ln);
        else if(se(bp,"t_phys")) t_phys=final_double(ln);
        else if(se(bp,"v_phys")) v_phys=final_double(ln);
        else if(se(bp,"g_phys")) g_phys=final_double(ln);
        else if(se(bp,"nu_phys")) nu_phys=final_double(ln);
        else if(se(bp,"c_L")) c_L=final_double(ln);
        else if(se(bp,"c_D")) c_D=final_double(ln);
        else if(se(bp,"box_size")) {
            lx=ly=lz=final_double(ln);
            if(lx<=0) fatal_error("Box size must be positive",1);
        } else if(se(bp,"box_size_xyz")) {
            lx=atof(next_token(ln));
            ly=atof(next_token(ln));
            lz=final_double(ln);
            if(lx<=0||ly<=0||lz<=0) fatal_error("Box sizes must be positive",1);
        } else if(se(bp,"wind_rms")) w_rms=final_double(ln);
        else if(se(bp,"wind_rms_phys")) w_rms_phys=final_double(ln);
        else if(se(bp,"wind_C")) {
            bp=next_token(ln);
            w_Cinv=se(bp,"frozen")?0:1/atof(bp);
            check_no_more(ln);
        } else if(se(bp,"wind_C_phys")) {
            bp=next_token(ln);
            if(se(bp,"frozen")) w_Cinv=0;
            else w_Cinv_phys=1/atof(bp);
            check_no_more(ln);
        } else if(se(bp,"wind_alpha")) w_alpha=final_double(ln);
        else if(se(bp,"gpr_rs")) gpr_rs=final_int(ln);
        else if(se(bp,"gpr_ts")) gpr_ts=final_int(ln);
        else if(se(bp,"gpr_dist_pad")) gpr_dist_pad=final_double(ln);
        else if(se(bp,"gpr_ker_tol")) {
            gpr_ker_tol=final_double(ln);
            if(gpr_ker_tol<0) fatal_error("The gpr_ker_tol should be positive",1);
        } else if(se(bp,"gpr_k_rt_param")) {
            gpr_rs=atoi(next_token(ln));
            gpr_ts=atoi(next_token(ln));
            gpr_dist_pad=final_double(ln);
        } else if(se(bp,"gpr_k_r_param")) {
            gpr_rs=atoi(next_token(ln));
            gpr_dist_pad=final_double(ln);
        } else if(se(bp,"glider_ts_pad")) gl_pad=final_double(ln);
        else if(se(bp,"wind_ts_pad")) tf_pad=final_double(ln);
        else if(se(bp,"ctl_duration")) c_dur=final_double(ln);
        else if(se(bp,"ctl_duration_phys")) c_dur_phys=final_double(ln);
        else if(se(bp,"bank_angle_step")) bstep=final_double(ln);
        else if(se(bp,"duration")) duration=final_double(ln);
        else if(se(bp,"duration_phys")) duration_phys=final_double(ln);
        else if(se(bp,"wind_model")) {
            bp=next_token(ln);
            if(se(bp,"gpr")) wmodel=wm_gpr;
            else if(se(bp,"full_linear")) wmodel=wm_full_linear;
            else if(se(bp,"full_cubic")) wmodel=wm_full_cubic;
            else fatal_error("Wind model type not understood\n",1);
            check_no_more(ln);
        } else if(se(bp,"integration_type")) {
            bp=next_token(ln);
            if(se(bp,"euler")) itype=it_euler;
            else if(se(bp,"improv_e")) itype=it_improv_e;
            else fatal_error("Integration type not understood\n",1);
            check_no_more(ln);
        } else if(se(bp,"planning")) {
            bp=next_token(ln);
            if(se(bp,"zero")) ptype=pt_zero;
            else if(se(bp,"random")) ptype=pt_random;
            else if(se(bp,"mcts")) ptype=pt_mcts;
            else fatal_error("Planning type not understood\n",1);
            check_no_more(ln);
        } else if(se(bp,"snapshots")) {
            bp=next_token(ln);
            while(bp!=NULL) {
                if(se(bp,"glider_xyz")) fflags|=16;
                else if(se(bp,"glider_mb")) fflags|=32;
                else if(se(bp,"glider_lb")) fflags|=64;
                else if(se(bp,"wind_modes")) fflags|=128;
                else fatal_error("Snapshot type not understood",1);
                bp=strtok(NULL," \t\n");
            }
        } else if(se(bp,"wind_slice")) {
            fflags|=8;
            int va=0,va2;

            // Read the direction of the cross-section to output
            bp=next_token(ln);
            if(se(bp,"y")) va=1;
            else if(se(bp,"z")) va=2;
            else if(!se(bp,"x"))
                fatal_error("First argument of wind_slice should be x, y, or z",1);

            // Read the wind component to output
            bp=next_token(ln);
            if(se(bp,"wy")) va|=4;
            else if(se(bp,"wz")) va|=8;
            else if(!se(bp,"wx"))
                fatal_error("Second argument of wind_slice should be wx, wy, or wz",2);

            // Read the grid index of the cross-section
            va2=final_int(ln);
            if(va2<0||va2>=(1>>28)) fatal_error("Index of wind_slice out of range",3);
            w_cs.push_back((va2<<4)|va);
        } else if(se(bp,"summary")) {
            bp=next_token(ln);
            while(bp!=NULL) {
                if(se(bp,"climb_stats")) fflags|=1;
                else if(se(bp,"en_components")) fflags|=2;
                else if(se(bp,"checksums")) fflags|=4;
                else fatal_error("Summary type not understood",1);
                bp=strtok(NULL," \t\n");
            }
        } else if(se(bp,"wind_cor")) {
            fflags|=256;
            wic_tcut=final_double(ln);
        } else if(se(bp,"wind_cor_phys")) {
            fflags|=256;
            wic_tcut_phys=final_double(ln);
        } else if(se(bp,"mcts_info")) {
            fflags|=512;
            mti_tcut=final_double(ln);
        } else if(se(bp,"mcts_info_phys")) {
            fflags|=512;
            mti_tcut_phys=final_double(ln);
        } else if(se(bp,"path_interval")) {
            fflags|=1024;
            path_interval=final_int(ln);
            if(path_interval<=0) fatal_error("Invalid value of path_interval",1);
        } else {
            fprintf(stderr,"Keyword '%s' not understood on line %d\n",bp,ln);
            exit(1);
        }
        ln++;
    }
    delete [] buf;
    fclose(f);

    // Check that all basic integer class parameters have been specified
    if(num_trials==-1) fatal_error("Number of trials not set",1);
    if(gpt==-1) fatal_error("Gliders per trial not set",1);
    if(nx==-1||ny==-1||nz==-1) fatal_error("Number of Fourier modes not set",1);
    if(br_min==1&&br_max==0) fatal_error("Bank angle range not set",1);
    if(bc_min==1&&bc_max==0) fatal_error("Bank angle control range not set",1);

    // Check that all basic float class parameters have been specified
    check_invalid(nu_phys,"nu_phys");
    check_invalid(c_L,"c_L");
    check_invalid(c_D,"c_D");
    check_invalid(lx,"lx");
    check_invalid(ly,"ly");
    check_invalid(lz,"lz");
    check_invalid(bstep,"bank_angle_step");

    // Check physical scales, and fill in others
    int nscales=(l_phys<0?1:0)+(t_phys<0?1:0)+(v_phys<0?1:0)+(g_phys<0?1:0);
    if(nscales!=2) fatal_error("Precisely two physical scales must be given",1);
    if(l_phys<0) {
        l_phys=t_phys<0?v_phys*v_phys/g_phys:
              (v_phys<0?v_phys*t_phys:g_phys*t_phys*t_phys);
    }
    if(t_phys<0) {
        t_phys=v_phys<0?sqrt(l_phys/g_phys):l_phys/v_phys;
    }
    if(v_phys<0) v_phys=l_phys/t_phys;
    if(g_phys<0) g_phys=v_phys/t_phys;

    // Convert any physical values into actual values
    if(duration_phys>=0) duration=duration_phys/t_phys;
    check_invalid(duration,"duration");
    if(c_dur_phys>=0) c_dur=c_dur_phys/t_phys;
    check_invalid(c_dur,"c_dur");
    if(w_Cinv_phys>=0) w_Cinv=w_Cinv_phys*t_phys;
    if(w_rms_phys>=0) w_rms=w_rms_phys/v_phys;

    // Check that any fluid cross-sections are within range
    for(unsigned int o=0;o<w_cs.size();o++) {
        int ind=w_cs[o]>>4,dir=w_cs[o]&3;
        if(ind>=(dir==0?nx:(dir==1?ny:nz)))
            fatal_error("Fluid cross-section index larger than box size",1);
    }

    // Set up the glider model
    gm=new glider_model(c_L,c_D,br_min,br_max,bc_min,bc_max,bstep);

    // Check wind parameters
    if(w_alpha<0) {
        if(w_rms<0) fatal_error("Need to specify either w_alpha or w_rms",1);
        calculate_wind_param(false);
    } else {
        if(w_rms>=0) fatal_error("Only one of w_alpha or w_rms should be specified",1);
        calculate_wind_param(true);
    }
    check_invalid(w_Cinv,"w_Cinv");

    // Check MCTS, and compute the z velocity scale for normalizing MCTS
    // scores
    if(ptype==pt_unset) fatal_error("No prediction model given",1);
    if(ptype==pt_mcts&&(n_mcts==-1||mcts_depth==-1||mcts_ex_fac<0))
        fatal_error("MCTS parameters not fully specified",1);

    // Check GPR parameters
    if(wmodel==wm_unset) fatal_error("No wind model given",1);
    if(wmodel==wm_gpr) {
        check_invalid(gpr_dist_pad,"gpr_dist_pad");
        if(gprm==-1||gs_pc==-1||gpr_rs==-1||(!frozen()&&gpr_ts==-1))
            fatal_error("GPR parameters not fully specified",1);
    } else {
        if(fflags&4) fatal_error("Checksums only valid for GPR prediction",1);
        if(!frozen()&&h_segs==-1)
            fatal_error("Hermite segments must be specified",1);
    }

    // Set up the wind correlation computation
    if(fflags&256) {
        if(wmodel!=wm_gpr) fatal_error("Wind correlation computation only works for GPR model",1);
        if(ptype!=pt_mcts) fatal_error("Wind correlation computation only works for MCTS",1);
        if(!frozen()) fatal_error("Wind correlation computation only valid for frozen field",1);
        if(wic_tcut_phys>=0) wic_tcut=wic_tcut_phys/t_phys;
    }

    // Check the MCTS tree info is consistent
    if(fflags&512) {
        if(ptype!=pt_mcts) fatal_error("Tree info requires MCTS option",1);
        if(mti_tcut_phys>=0) mti_tcut=mti_tcut_phys/t_phys;
    }

    // Check that MCTS is enabled for path output
    if(fflags&1024&&ptype!=pt_mcts) fatal_error("Path output requires MCTS option",1);

    // Check integration model
    if(itype==it_unset) fatal_error("No integration model specified",1);
}

/** Chooses the simulation timestep based on the timestep restrictions
 * from the turbulent fluid simulation and the glider model. Computes
 * other constants needed for integration.
 * \param[in] verbose whether to output information about the calculations. */
void fileinfo::select_timestep(bool verbose) {
    const double fmax=std::numeric_limits<double>::max(),
                 oneminus=(1-16*std::numeric_limits<double>::epsilon());
    bool tf_used;

    // Compute the number of integration intervals to divide a control event
    // into, based on dividing by the least common multiple of the memory
    // storage frequency and the output frequency
    i_pc=wmodel!=wm_gpr?out_pc:lcm(gs_pc,out_pc);
    i_dur=c_dur/i_pc;
    out_dur=c_dur/out_pc;

    // Compute the timestep restrictions for the turbulent fluid model and the
    // glider model. Multiply each one by a padding factor.
    double tf_max=frozen()?fmax:max_tf_timestep(),
           tf_ts=frozen()?fmax:tf_max*tf_pad,
           gl_max=gm->max_timestep(),gl_ts=gl_max*gl_pad,dt_orig;

    // Choose the minimum of the two padded timesteps
    if(tf_ts<gl_ts) {tf_used=true;dt_orig=tf_ts;}
    else {tf_used=false;dt_orig=gl_ts;}

    // Compute adjusted timestep so that an integer number will cover the
    // integration interval. When the glider correlation computation is
    // enabled, adjust the timestep to be a multiple of 10.
    int pre=fflags&256?10/gcd(10,i_pc):1;
    steps=pre*int((1./pre)*oneminus*i_dur/dt_orig+1);
    dt=i_dur/steps;

    // Compute the total number of control events, integration intervals,
    // and timesteps
    int c_events=int(oneminus*duration/c_dur+1);
    snaps=c_events*out_pc;
    intervals=c_events*i_pc;

    // Compute the frequencies for GPR sampling and output
    gs_freq=i_pc/gs_pc;
    out_freq=i_pc/out_pc;

    // If requested, print information about these calculations
    if(verbose) {
        double T_phys=t_phys*1e3;

        // Memory measurements are only relevant when GPR is in used
        if(wmodel!=wm_gpr) puts("# Memory measurements per control event : [not used]");
        else printf("# Memory measurements per control event : %d\n",gs_pc);

        // Print the rest of the information
        printf("# Outputs per control event             : %d\n"
               "# Substeps per control event            : %d\n"
               "# Control event duration                : %g (%g s)\n"
               "# Integration interval duration         : %g (%g s)\n#\n",
               out_pc,i_pc,c_dur,c_dur*t_phys,i_dur,i_dur*t_phys);
        if(frozen())
            puts("# Max fluid timestep                    : [frozen]\n"
                 "# Fluid timestep padding factor         : [frozen]\n"
                 "# Fluid timestep                        : [frozen]\n#");
        else
            printf("# Max fluid timestep                    : %g (%g ms)\n"
                   "# Fluid timestep padding factor         : %g\n"
                   "# Fluid timestep                        : %g (%g ms)%s\n#\n",
                   tf_max,tf_max*T_phys,tf_pad,tf_ts,tf_ts*T_phys,tf_used?" <-- use this":"");
        printf("# Max glider timestep                   : %g (%g ms)\n"
               "# Glider timestep padding factor        : %g\n"
               "# Glider timestep                       : %g (%g ms)%s\n#\n"
               "# Adjusted timestep                     : %g (%g ms)\n"
               "# Timesteps per integration interval    : %d\n#\n"
               "# Simulation duration                   : %g (%g s)\n"
               "# Total control events                  : %d\n"
               "# Total output snapshots                : %d [+1]\n"
               "# Total integration intervals           : %d\n"
               "# Total timesteps                       : %d\n",
               gl_max,gl_max*T_phys,gl_pad,gl_ts,gl_ts*T_phys,
               tf_used?"":" <-- use this",dt,dt*T_phys,steps,duration,
               duration*t_phys,c_events,snaps,intervals,steps*intervals);
    }
}

/** Prints information about all of the constants stored within the class.
 * \param[in] fp a file handle to write to. */
void fileinfo::print_info(FILE *fp) {

    // Output information about physical scales, the overall simulation setup,
    // and the glider parameters
    double gspeed=1/sqrt(gm->s);
    fprintf(fp,"# Physical scales\n"
               "# -----------------------\n"
               "# Length scale          : %g m\n"
               "# Gravity scale         : %g m/s^2\n"
               "# Time scale            : %g s\n"
               "# Velocity scale        : %g m/s\n#\n"
               "# Simulation setup\n"
               "# -----------------------\n"
               "# Number of trials      : %d\n"
               "# Gliders per trial     : %d\n"
               "# Duration              : %g (%g s)\n"
               "# Random seed           : %lu\n#\n"
               "# Glider information\n"
               "# -----------------------\n"
               "# Lift coefficient      : %g\n"
               "# Drag coefficient      : %g\n"
               "# Glider speed          : %g (%g m/s)\n"
               "# Control duration      : %g (%g s)\n"
               "# Bank angles (deg)     : ",
               l_phys,g_phys,t_phys,v_phys,
               num_trials,gpt,duration,duration*t_phys,base_seed,
               gm->c_L,gm->c_D,gspeed,gspeed*v_phys,
               c_dur,c_dur*t_phys);

    // Print the bank angle information, showing a range of the possible values
    print_range(fp,gm->br_min,gm->br_max);
    fputs("\n# Bank changes (deg)    : ",fp);
    print_range(fp,gm->bc_min,gm->bc_max);

    // Print information about the turbulent wind field
    fprintf(fp,"\n# Integration method    : %s\n#\n"
               "# Wind field information\n"
               "# -----------------------\n"
               "# Modes                 : %d x %d x %d\n"
               "# Box dimensions        : %g x %g x %g\n"
               "# Box dimensions (SI)   : (%g m) x (%g m) x (%g m)\n"
               "# alpha parameter       : %g\n"
               "# Wind RMS              : %g (%g m/s)\n"
               "# C parameter           : ",
               s_integration_type(),nx,ny,nz,lx,ly,lz,
               lx*l_phys,ly*l_phys,lz*l_phys,w_alpha,w_rms,w_rms*v_phys);

    // Print information about the alpha parameter if it is in use
    if(frozen()) fputs("[frozen]",fp);
    else fprintf(fp,"%g (%g s)",1/w_Cinv,t_phys/w_Cinv);

    // Print information about the wind prediction model
    fprintf(fp,"\n#\n# Wind prediction model\n"
               "# -----------------------\n"
               "# Model type            : %s\n",
               s_wind_model());

    // Print extra parameters needed for the selected prediction model
    if(wmodel==wm_gpr) {
        fprintf(fp,"# GPR memory            : %d\n"
                   "# Samples per ctl. dur. : %d\n",gprm,gs_pc);
        frozen()?fprintf(fp,"# Kernel sample grid    : %d\n",gpr_rs)
                :fprintf(fp,"# Kernel sample grid    : %d by %d\n",gpr_rs,gpr_ts);
        fprintf(fp,"# Distance pad factor   : %g\n"
                   "# Kernel update tol.    : %g\n",gpr_dist_pad,gpr_ker_tol);
    } else if(!frozen()) {
        fprintf(fp,"# Hermite segments      : %d\n",h_segs);
    }

    // Print the planning model type
    fprintf(fp,"#\n# Planning model\n"
            "# -----------------------\n"
            "# Model type            : %s\n",
            s_planning_type());

    // Print extra information about the MCTS, if it is in use
    if(ptype==pt_mcts) {
        fprintf(fp,"# MCTS trials           : %d\n"
                "# MCTS tree depth       : %d\n"
                "# MCTS exploration fac. : %g\n",
                n_mcts,mcts_depth,mcts_ex_fac);
    }

    // Print information about the output types
    fprintf(fp,"#\n# Output information\n"
            "# -----------------------\n"
            "# Outputs per ctl. dur. : %d\n"
            "# Snapshot types        :",out_pc);
    if((fflags&240)==0) fputs(" [none]",fp);
    else {
        if(fflags&16) fputs(" glider_xyz",fp);
        if(fflags&32) fputs(" glider_lb",fp);
        if(fflags&64) fputs(" glider_fb",fp);
        if(fflags&128) fputs(" wind_modes",fp);
    }

    // Print wind cross-sections
    fputs("\n# Wind cross-sections   :",fp);
    if(fflags&8) {
        for(unsigned int o=0;o<w_cs.size();o++) {
            int ind=w_cs[o]>>4,dir=w_cs[o]&3,cmp=(w_cs[o]>>2)&3;
            fprintf(fp," (w%c,%c,%d)",'x'+cmp,'x'+dir,ind);
        }
    } else fputs(" [none]",fp);

    // Print summary types
    fputs("\n# Summary types         :",fp);
    if((fflags&7)==0) fputs(" [none]\n",fp);
    else {
        if(fflags&1) fputs(" climb_stats",fp);
        if(fflags&2) fputs(" en_components",fp);
        if(fflags&4) fputs(" checksums",fp);
        fputc('\n',fp);
    }

    // Print information about the wind correlation computation
    if(fflags&256)
        fprintf(fp,"# Wind correl. cutoff   : %g (%g s)\n",wic_tcut,wic_tcut*t_phys);

    // Print information about the MCTS tree info
    if(fflags&512)
        fprintf(fp,"# MCTS tree info cutoff : %g (%g s)\n",mti_tcut,mti_tcut*t_phys);
}

/** Prints a range of glider angles.
 * \param[in] fp the file handle to write to.
 * \param[in] (rmin,rmax) the integer range for the angles, which should be
 *                        multipled by the glider banking angle step size. */
void fileinfo::print_range(FILE *fp,int rmin,int rmax) {
    for(int i=rmin;i<rmax;i++) fprintf(fp,"%g,",i*gm->bstep);
    fprintf(fp,"%g",rmax*gm->bstep);
}

/** Checks that there are no subsequent values.
 * \param[in] ln the current line number. */
void fileinfo::check_no_more(int ln) {
    if(strtok(NULL," \t\n")!=NULL) {
        fprintf(stderr,"Too many arguments at input line %d\n",ln);
        exit(1);
    }
}

/** Finds the next token in a string, and if none is availble, gives an error
 * message.
 * \param[in] ln the current line number. */
char* fileinfo::next_token(int ln) {
    char *temp=strtok(NULL," \t\n");
    if(temp==NULL) {
        fprintf(stderr,"Not enough arguments at input line %d\n",ln);
        exit(1);
    }
    return temp;
}

/** Checks that a parameter is a valid positive value.
 * \param[in] val the value to check.
 * \param[in] p the name of the value. */
void fileinfo::check_invalid(double val,const char *p) {
    if(val<0) {
        fprintf(stderr,"Value of %s either invalid or not set\n",p);
        exit(1);
    }
}

/** Sets either the wind alpha parameter or the wind RMS in order to be consistent
 * with the other. */
inline void fileinfo::calculate_wind_param(bool set_rms) {
    double s=0,facx=2*M_PI/lx,facy=2*M_PI/ly,facz=2*M_PI/lz;

#pragma omp parallel for reduction(+:s)
    for(int k=0;k<nz;k++) {
        double w=sqr(facz*(k>nz/2?nz-k:k));
        for(int j=0;j<ny;j++) {
            double ww=w+sqr(facy*(j>ny/2?ny-j:j));
            for(int i=0;i<((nx>>1)+1);i++) if(f_mode(i,j,k)!=0)
                s+=pow(ww+sqr(facx*i),-11/6.);
        }
    }

    // Fill in either the wind RMS or the alpha parameter
    if(set_rms) w_rms=sqrt(0.25/M_PI*facx*facy*facz*s*w_alpha);
    else w_alpha=4*M_PI*w_rms*w_rms/(facx*facy*facz*s);
}

/** Computes the maximum integration timestep for the turbulent fluid
 * simulation, based on the relaxation timescale of the fastest mode.
 * \return The calculated timestep. */
inline double fileinfo::max_tf_timestep() {
    double fx=nx/lx,fy=ny/ly,fz=nz/lz;
    return 2*pow(M_PI,-2/3)/w_Cinv*pow(fx*fx+fy*fy+fz*fz,-1/3.);
}

/** Computes the greatest common divisor of two integers.
 * \param[in] (a,b) the two integers.
 * \return The greatest common divisor. */
int fileinfo::gcd(int a,int b) {
    return b==0?a:gcd(b,a%b);
}
