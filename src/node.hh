#ifndef NODE_HH
#define NODE_HH

struct node {
    /** The number of visits. */
    int n;
    /** The index of the child node. If set to zero, no child has been created.
     */
    int c;
    /** The index of the move that this node represents. */
    short act;
    /** The number of branches at the child node, only used if the child node
     * has been created. */
    short nact;
    /** The expected win. */
    float w;
    inline void init(short act_) {
        n=0;c=0;act=act_;w=0.;
    }
    inline float weight() {
        return w/n;
    }
    inline void update(float w_) {
        n++;
        w+=w_;
    }
};

#endif
