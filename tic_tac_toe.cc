#include "tic_tac_toe.hh"

#include <cstdio>

tic_tac_toe::tic_tac_toe() {
    reset();
}

void tic_tac_toe::print() {
    printf("%s|%s|%s\n"
           "-+-+-\n"
           "%s|%s|%s\n"
           "-+-+-\n"
           "%s|%s|%s\n\n",
           c(0),c(1),c(2),c(3),c(4),
           c(5),c(6),c(7),c(8));
}

inline const char* tic_tac_toe::c(int i) {
    return board[i]==0?" ":(board[i]==1?"o":"x");
}

/** A move on the board.
 * \param[in] k grid index on the board.
 * \param[in] p player id (1 or 2).
 * \return the game status: win(+1), draw(0.5), or continue(0). */
float tic_tac_toe::play(int k,int p) {
    // 678
    // 345
    // 012

    // need to check if the move is valid??
    board[k]=p;
    if(full()) return 0.5;

    // check column
    int j=k%3;
    if(board[j]==p && board[j+3]==p && board[j+6]==p) return 1;

    // check row
    int i=k-j;
    if(board[i]==p && board[i+1]==p && board[i+2]==p) return 1;

    // Check diag
    if(board[4]==p){
        if(board[2]==p && board[6]==p) return 1;
        if(board[0]==p && board[8]==p) return 1;
    }

    return 0;

    // Return true if won
    // Return false otherwise
}

bool tic_tac_toe::full() {
    int s=0;
    for(int k=0;k<9;k++) if(board[k]!=0) s++;
    return s==9;
}
