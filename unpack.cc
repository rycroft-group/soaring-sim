#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>

#include "common.hh"

/** A maximum safe limit on the filename size. */
const int filename_max_size=8192;

/** The maximum number of gliders per trial to allow, to prevent very large
 * array allocations. */
const int gpt_max=131072;

inline void syntax_message() {
    fputs("./unpack [-w] <output_dir> <trial> <glider> <output_file>\n\n"
          "The -w option outputs the wind velocity along the glider's\n"
          "trajectory. It requires that the larger binary file was saved\n"
          "during the simulation.\n",stderr);
}

/** Checks whether a given file exists.
 * \param[in] fname the name of the file to check.
 * \return True if the file exists, false otherwise. */
bool file_exists(char *fname) {
    struct stat s;
    return stat(fname,&s)==0;
}

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc<5||argc>6) {
        syntax_message();
        return 1;
    }

    // Check for the wind flag
    bool wind;
    if(argc==6) {
        if(strcmp(argv[1],"-w")!=0) {
            syntax_message();
            return 1;
        }
        wind=true;
    } else wind=false;

    // Allocate a temporary array for assembling the input filename
    int len=strlen(argv[argc-4]),trial=atoi(argv[argc-3]),gnum=atoi(argv[argc-2]),gpt,k=0;
    if(len>filename_max_size) fatal_error("Filename too long",1);
    char *tfile=new char[len+128];

    // Check for a binary file that has the required information
    int b;
    if(wind) {

        // If wind information is requested, then only the larger binary file
        // can be used
        sprintf(tfile,"%s/glider_lb.%d",argv[argc-4],trial);
        if(!file_exists(tfile))
            fatal_error("Large binary file with wind data not available",1);
        b=6;
    } else {

        // If wind information isn't needed, then either binary file can be
        // used, but the minimal one is preferable because it is smaller
        sprintf(tfile,"%s/glider_mb.%d",argv[argc-4],trial);
        if(file_exists(tfile)) b=3;
        else {
            sprintf(tfile,"%s/glider_lb.%d",argv[argc-4],trial);
            if(!file_exists(tfile))
                fatal_error("Neither binary file available",1);
            b=6;
        }
    }
    FILE *fp=safe_fopen(tfile,"rb"),*outf;

    // Read the header information, with the number of gliders per file and the
    // output duration
    double out_dur;
    if((fread(&gpt,sizeof(int),1,fp)!=1)
     ||(fread(&out_dur,sizeof(double),1,fp)!=1))
        fatal_error("Error reading file header",1);
    if(gpt>=gpt_max) fatal_error("Gliders per trial exceeds a safe maximum",1);
    if(gnum<0||gnum>=gpt) fatal_error("Glider number out of range",1);

    // Set up temporary memory and pointers
    size_t rec=b*sizeof(float)+sizeof(short),recf=rec*gpt;
    char *ibuf=new char[recf];
    float *f=reinterpret_cast<float*>(ibuf+gnum*rec);
    short *pbank=reinterpret_cast<short*>(f+b);

    // Output the unpacked data
    bool std=strcmp(argv[argc-1],"-")==0;
    outf=std?stdout:safe_fopen(argv[argc-1],"w");
    while(fread(ibuf,recf,1,fp)==1) {
        wind?fprintf(outf,"%d %g %g %g %g %g %g %g %hd\n",
                     k,out_dur*k,*f,f[1],f[2],f[3],f[4],f[5],*pbank)
            :fprintf(outf,"%d %g %g %g %g %hd\n",k,out_dur*k,*f,f[1],f[2],*pbank);
        k++;
    }

    // Close any open files, and delete temporary memory
    fclose(fp);
    if(!std) fclose(outf);
    delete [] tfile;
    delete [] ibuf;
}
