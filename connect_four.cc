#include "connect_four.hh"
#include <cstdio>

connect_four_orig::connect_four_orig() {
    reset();
}

void connect_four_orig::print()
{
    printf(
	   "%s|%s|%s|%s|%s|%s|%s\n"
           "-+-+-+-+-+-+-\n"
           "%s|%s|%s|%s|%s|%s|%s\n"
           "-+-+-+-+-+-+-\n"
           "%s|%s|%s|%s|%s|%s|%s\n"
           "-+-+-+-+-+-+-\n"
           "%s|%s|%s|%s|%s|%s|%s\n"
           "-+-+-+-+-+-+-\n"
           "%s|%s|%s|%s|%s|%s|%s\n"
           "-+-+-+-+-+-+-\n"
           "%s|%s|%s|%s|%s|%s|%s\n"
	   "%d|%d|%d|%d|%d|%d|%d\n",
           c(35),c(36),c(37),c(38),c(39),c(40),c(41),
           c(28),c(29),c(30),c(31),c(32),c(33),c(34),
           c(21),c(22),c(23),c(24),c(25),c(26),c(27),
           c(14),c(15),c(16),c(17),c(18),c(19),c(20),
           c(7),c(8),c(9),c(10),c(11),c(12),c(13),
           c(0),c(1),c(2),c(3),c(4),c(5),c(6),
	   0,1,2,3,4,5,6);

}

inline const char* connect_four_orig::c(int i) {
    return board[i]==0?" ":(board[i]==1?"o":"x");
}

/** A move on the board.
 * \param[in] k grid index on the board.
 * \param[in] p player id (1 or 2).
 * \return the game status: win(+1), draw(0.5), or continue(0). */
float connect_four_orig::play(int k,int p) {

    // board

    // 35 36 37 38 39 40 41
    // ... ...
    // 14 15 16 17 18 19 20
    // 7  8  9  10 11 12 13
    // 0  1  2  3  4  5  6

    /// Play
    int i=0;
    while(i<6){
        int pos=i*7+k;
        if(board[pos]==0){
            board[pos]=p;
	  //  printf("played i is %d, k is %d\n",i,k);
            //rest_move--;
            break;
        }
        i++;
    }
//printf("you play ith row %d, kth col %d \n",i,k);
    // Have a new move in ith row, kth col
    // Check game status and return rewards
   // if(full()) return 0.5;

    // Check horizontal
    int count=0;
    int row,col;
    for(col=0;col<7;col++){
        if(board[i*7+col]==p) count++;
        else count=0;
//	printf("the col is %d, the count is %d\n",col,count);
	if(count==4) return 1;
    }

    // Check vertical
    count=0;
    for(row=0;row<6;row++){
        if(board[row*7+k]==p) count++;
        else count=0;
	if(count==4) return 1;

    }

    // Check diagonal
    count=0;
    if(i<k){
        row=0;
        col=k-i;
    }
    else{
        row=i-k;
        col=0;
    }

   while(row<6&&col<7){
        if(board[row*7+col]==p) {count++;}
        else {count=0;}
	if(count==4) return 1;
        row++;col++;

    }

    count=0;
    if(i<6-k){
        row=0;
        col=k+i;
    }
    else{
        row=i-6+k;
        col=6;
    }
    while(row<6&&col>-1){
        if(board[row*7+col]==p) {count++;}
        else {count=0;}
	// printf("the col is %d,row is %d,the count is %d\n",col,row,count);
	if(count==4) return 1;
	row++;col--;

    }
    return 0;
}

bool connect_four_orig::full() {
    int s=0;
    for(int k=0;k<42;k++) if(board[k]!=0) s++;
    return s==42;
}

/** Check if the move is valid. */
bool connect_four_orig::valid_move(int k){
  return board[35+k]==0;
}

/** Remove a step. */
void connect_four_orig::remove(int k){
  int i=5;
  while(i>-1){
    if(board[i*7+k]!=0){
      board[7*i+k]=0;
      break;
    }
    i--;

  }

}
