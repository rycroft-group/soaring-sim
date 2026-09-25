#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "common.hh"
#include "fileinfo.hh"

/** A maximum safe limit on the filename size. */
const int filename_max_size=8192;

/** The maximum number of gliders per trial to allow, to prevent very large
 * array allocations. */
const int gpt_max=131072;

/** A small multiple of machine epsilon used as a tolerance in arithmetic. */
const double tol=16*std::numeric_limits<double>::epsilon();

/** Tolerances to ensure correct rounding in computing the output frame range
 * to use. */
const double oneminus=1-tol,oneplus=1+tol;

/** Checks if two floating point numbers differ significantly. */
inline bool far(double a,double b) {
    return fabs(a-b)>tol*(fabs(a)+fabs(b));
}

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc<3||argc>4) {
        fputs("./climb_rate <input_file> <t_lo> [<t_hi>]\n\n"
              "Calculates the climb rate of individual gliders based on the trajectory\n"
              "data over a time interval t_lo<=t<=t_hi. If t_hi is omitted, the"
              "calculation is performed for t>=t_hi",stderr);
        return 1;
    }

    // Load the configuration file. This is needed to obtain information about
    // the time scales in the simulation.
    fileinfo fi(argv[1]);

    // Compute the first output frame to use
    double t_lo=atof(argv[2]),t_hi,fout_dur,
           out_dur=fi.c_dur/fi.out_pc,odp=out_dur*fi.t_phys;
    int snaps=int(oneminus*fi.duration/fi.c_dur+1)*fi.out_pc,
        k_lo=int(oneminus*t_lo/odp+1),k_hi;
    if(k_lo>snaps) fatal_error("t_min out of range",1);

    // Compute the last output frame to use
    if(argc==4) {
        t_hi=atof(argv[3]);
        k_hi=int(oneplus*t_hi/odp);
        if(k_hi>snaps) fatal_error("t_max out of range",1);
    } else k_hi=snaps;

    // Allocate memory for computing the regression coefficients
    int j,k,l,n=k_hi-k_lo+1,gpt=fi.gpt,fgpt;
    double *s=new double[2*gpt],
           fac=1./(n*(n+1)),fac2=fac/out_dur,
           o=(4*n-2)*fac,p=-6*fac,q=-6*fac2,r=12./(n-1)*fac2,
           nor=1./(out_dur*(n-1));
    float *v=new float[gpt];

    // Allocate the buffer for reading the file
    size_t rec=3*sizeof(float)+sizeof(short),recf=rec*gpt;
    char *buf=new char[strlen(fi.filename)+128],*ibuf=new char[recf];

    // Open the output file and output a header line about the range of output
    // frames that are used
    sprintf(buf,"%s/cstats_indiv",fi.filename);
    FILE *outf=safe_fopen(buf,"w");
    if(argc==4) fprintf(outf,"Computed for %g s <= t <= %g s (Frames %d to %d)\n",t_lo,t_hi,k_lo,k_hi);
    else fprintf(outf,"# Computed for t >= %g s (Frames %d to %d)\n",t_lo,k_lo,k_hi);

    // Loop over the different trajectory files
    for(j=0;j<fi.num_trials;j++) {

        // Open the glider trajectory file and check for header consistency
        sprintf(buf,"%s/glider_mb.%d",fi.filename,j);
        FILE *fp=safe_fopen(buf,"rb");
        if((fread(&fgpt,sizeof(int),1,fp)!=1)
         ||(fread(&fout_dur,sizeof(double),1,fp)!=1))
                fatal_error("Error reading file header",1);
        if(fgpt!=gpt||far(fout_dur,out_dur)) fatal_error("Parameter mismatch in header",1);

        // Throw out the first section of the file
        if(fseek(fp,recf*k_lo,SEEK_CUR)!=0)
            fatal_error("Can't find start of snapshot data to use",1);

        // Loop over the snapshots to read
        for(int i=0;i<2*gpt;i++) s[i]=0;
        for(k=0;k<n;k++) {
            if(fread(ibuf,recf,1,fp)!=1) {
                fprintf(stderr,"Error reading snapshot %d in file %d\n",k,j);
                exit(1);
            }

            // Compute each glider's contribution to the linear regression
            for(l=0;l<gpt;l++) {
                float *f=reinterpret_cast<float*>(ibuf+rec*l);
                s[2*l]+=f[2];
                s[2*l+1]+=k*f[2];
                if(k==0) v[l]=f[2];
            }
        }

        // Output the regression coefficients to the output file
        for(l=0;l<gpt;l++) {
            float *f=reinterpret_cast<float*>(ibuf+rec*l);
            fprintf(outf,"%d %d %g %g %g\n",j,l,
                    o*s[2*l]+p*s[2*l+1],q*s[2*l]+r*s[2*l+1],
                    nor*(double(f[2])-v[l]));
        }
        fclose(fp);
    }

    // Close output file and delete dynamically allocated memory
    delete [] ibuf;
    delete [] buf;
    delete [] v;
    delete [] s;
    fclose(outf);
}
