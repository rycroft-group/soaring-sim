#ifndef PATH_STORE_HH
#define PATH_STORE_HH

#include <cstdio>

#include "glider.hh"

struct path_data {
    /** The glider x position. */
    float rx;
    /** The glider y position. */
    float ry;
    /** The glider z position. */
    float rz;
    /** The x component of the predicted wind at the glider's position. */
    float wx;
    /** The y component of the predicted wind at the glider's position. */
    float wy;
    /** The z component of the predicted wind at the glider's position. */
    float wz;
    /** The glider's banking position. */
    short bank;
    inline void set(float rx_,float ry_,float rz_,float wx_,float wy_,float wz_,short bank_) {
        rx=rx_;ry=ry_;rz=rz_;
        wx=wx_;wy=wy_;wz=wz_;
        bank=bank_;
    }
    inline void write_binary(FILE *fp) {
        fwrite(&rx,sizeof(float),6,fp);
        fwrite(&bank,sizeof(short),1,fp);
    }
};

struct path_info {
    /** The counter for storing information within this class. */
    int c;
    /** The size of the path data array. */
    int n;
    /** The array of path data. */
    path_data* const d;
    path_info(int n_) : c(0), n(n_), d(new path_data[n]) {}
    /** The class destructor frees the dynamically allocated memory. */
    ~path_info() {
        delete [] d;
    }
    inline void store(glider &g,double *w) {
        d[c++].set(g.rx,g.ry,g.rz,*w,w[1],w[2],g.bank);
    }
    inline void output_and_reset(FILE *fp) {
        for(int i=0;i<n;i++) d[i].write_binary(fp);
        c=0;
    }
};

#endif
