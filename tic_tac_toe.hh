#ifndef TIC_TAC_TOE_HH
#define TIC_TAC_TOE_HH

class tic_tac_toe {
    public:
        int board[9];
        tic_tac_toe();
        void print();
	    float play(int k,int p);
        inline const char* c(int i);
        inline void reset() {
            for(int i=0;i<9;i++) board[i]=0;
        }
        bool full();
};

#endif
