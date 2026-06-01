#ifndef CONNECT_FOUR
#define CONNECT_FOUR

class connect_four_orig {
    public:
        int board[42];
        connect_four_orig();
        void print();
	    float play(int k,int p);
        inline const char* c(int i);
        inline void reset() {
            for(int i=0;i<42;i++) board[i]=0;
        }
        bool full();
        bool valid_move(int k);
        void remove(int k);
};

#endif
