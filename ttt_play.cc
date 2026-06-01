#include "connect_four.hh"
#include "mcts.hh"
#include "tic_tac_toe.hh"
int main() {
//
	connect_four_orig ttt;
//	tic_tac_toe ttt;
	mcts_orig mcts(7,ttt);
	bool w;
	int action;

	// test case
	/**
	printf("full? %d\n",ttt.full());
	w=ttt.play(0,1);
	printf("w is %d\n",w);
	 w=ttt.play(1,1);
        printf("w is %d\n",w);
	 w=ttt.play(2,2);
        printf("w is %d\n",w);
	 w=ttt.play(3,1);
	 ttt.play(4,2);ttt.play(5,1);ttt.play(6,1);
	 ttt.play(0,2);ttt.play(1,2);ttt.play(2,1);ttt.play(3,1);ttt.play(4,1);ttt.play(5,2);ttt.play(6,1);
	  ttt.play(0,1);ttt.play(1,1);ttt.play(2,2);ttt.play(3,2);ttt.play(4,2);ttt.play(5,1);ttt.play(6,1);

	   ttt.play(0,2);ttt.play(1,2);ttt.play(2,1);ttt.play(3,1);ttt.play(4,2);ttt.play(5,1);ttt.play(6,2);

	    ttt.play(0,2);ttt.play(1,1);ttt.play(2,2);ttt.play(3,2);ttt.play(4,1);ttt.play(5,1);ttt.play(6,1);

	     ttt.play(0,2);ttt.play(1,2);ttt.play(2,2);w=ttt.play(3,1);ttt.play(4,1);ttt.play(5,2);ttt.play(6,2);

        printf("w is %d\n",w);
	ttt.print();
	ttt.reset();
	ttt.play(3,1);ttt.play(3,1);ttt.play(3,1);ttt.play(3,2);ttt.play(3,1);ttt.play(3,2);
	ttt.play(4,1);ttt.play(4,1);ttt.play(4,2);ttt.play(4,1);ttt.play(4,2);
	ttt.play(5,1);ttt.play(5,1);ttt.play(5,1);ttt.play(5,2);
	ttt.play(6,2);ttt.play(6,2);w=ttt.play(6,2);
        printf("w is %d\n",w);
	ttt.print();
*/

	while(true) {

	// Play computer move
	puts("Computer's move:");
		mcts.reset();
		for(int t=0;t<10000;t++) mcts.simulate_path(1,t);
		printf("the size is %d\n",mcts.s);
		action=mcts.pick_best(0,1);
		w=ttt.play(action,1);
        ttt.print();
	printf("Computer played %d\n",action);

        if(w==1) {
            puts("Sorry ... better luck next time");
            return 0;
        }
	if(ttt.full()){
	    puts("Aha! There is a draw!");
	    return 0;
	}

        // Ask the user to make a move, checking that it is in range
	puts("What is your move?");
//	 printf("Recommand move:  %d\n",mcts.pick_best(0,2));
        scanf("%d",&action);
        while(action<0||action>6) {
            puts("Choose a number between 0 and 6:\n\n0|1|2"
                 "\n-+-+-\n3|4|5\n-+-+-\n6|7|8\n\n");
            scanf("%d",&action);
        }
	while(ttt.board[35+action]!=0){
		puts("The column is full, please take another move");
		scanf("%d",&action);
	}
        // Play the user move
        w=ttt.play(action,2);
        ttt.print();
        if(w==1) {
            puts("Congratulations!!!");
            return 0;
        }
	if(ttt.full()){
            puts("Aha! There is a draw!");
	    return 0;
        }
	}

}
