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

/** Checks whether a given file exists.
 * \param[in] fname the name of the file to check.
 * \return True if the file exists, false otherwise. */
bool file_exists(char *fname) {
    struct stat s;
    return stat(fname,&s)==0;
}

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc!=7) {
        fputs("./get_paths <output_dir> <trial> <glider> <step> <path> <output_file>\n\n"
              "If <path> is set to '-', then all paths are outputted. If <output_file>\n"
              "is set to '-', then the program writes to standard output.\n",stderr);
        return 1;
    }

    // Allocate a temporary array for assembling the input filename
    int len=strlen(argv[1]),trial=atoi(argv[2]),gnum=atoi(argv[3]),
        step=atoi(argv[4]),gpt,n_mcts,mcts_depth,out_pc,jlo,jhi;
    if(step<=0) fatal_error("<step> argument must be greater than zero",1);
    if(len>filename_max_size) fatal_error("Filename too long",1);
    char *tfile=new char[len+128];

    // Open the path information file
    snprintf(tfile,len+128,"%s/path_info.%d",argv[1],trial);
    FILE *fp=safe_fopen(tfile,"rb"),*outf;

    // Read the header information, with the number of gliders per file and the
    // output duration
    double p_dur,out_dur;
    if((fread(&gpt,sizeof(int),1,fp)!=1)
     ||(fread(&n_mcts,sizeof(int),1,fp)!=1)
     ||(fread(&mcts_depth,sizeof(int),1,fp)!=1)
     ||(fread(&out_pc,sizeof(int),1,fp)!=1)
     ||(fread(&p_dur,sizeof(double),1,fp)!=1)
     ||(fread(&out_dur,sizeof(double),1,fp)!=1))
        fatal_error("Error reading file header",1);
    if(gpt>=gpt_max) fatal_error("Gliders per trial exceeds a safe maximum",1);
    if(gnum<0||gnum>=gpt) fatal_error("Glider number out of range",1);

    // Check for the path(s) to output
    if(strcmp(argv[5],"-")==0) {jlo=0;jhi=n_mcts;}
    else {
        jlo=atoi(argv[5]);
        if(jlo<0||jlo>=n_mcts) fatal_error("Path number out of range",1);
        jhi=jlo+1;
    }

    // Seek to the correct block of information in the file
    int q=1+mcts_depth*out_pc;
    size_t rec=6*sizeof(float)+sizeof(short),
           blo=rec*q*n_mcts;
    if(fseek(fp,blo*(gnum+gpt*(step-1)),SEEK_CUR)!=0)
        fatal_error("Error finding path data in file",1);

    // Read in the block of information
    char *ibuf=new char[blo];
    if(fread(ibuf,blo,1,fp)!=1)
        fatal_error("Error reading path data from file",1);

    // Output the path information to a file
    bool std=strcmp(argv[6],"-")==0;
    outf=std?stdout:safe_fopen(argv[6],"w");
    double base_dur=p_dur*step;
    for(int i=0;i<q;i++) {
        fprintf(outf,"%d %g",i,base_dur+out_dur*i);
        for(int j=jlo;j<jhi;j++) {
            float *f=reinterpret_cast<float*>(ibuf+rec*(j*q+i));
            short *pbank=reinterpret_cast<short*>(f+6);
            fprintf(outf," %g %g %g %g %g %g %hd",
                    *f,f[1],f[2],f[3],f[4],f[5],*pbank);
        }
        fputc('\n',outf);
    }

    // Close any open files, and delete temporary memory
    fclose(fp);
    if(!std) fclose(outf);
    delete [] tfile;
    delete [] ibuf;
}
