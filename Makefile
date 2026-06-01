# Load the common configuration file
include ../config.mk

iflags=`gsl-config --cflags` $(fftw_iflags)
lflags=`gsl-config --libs` $(fftw_lflags)

objs=turb_fluid.o tf_grid.o tf_grid_mr.o glider.o common.o lp_solve.o k_func.o gpr.o tic_tac_toe.o mcts.o connect_four.o game_mcts.o game.o file_output.o fileinfo.o soaring_sim.o soaring_sim_io.o wind_correl.o
src=$(patsubst %.o,%.cc,$(objs))
execs=soar t_correl gl_conv gl_test t_vel_stats t_correl kf_calc gpr_test2 lapack_test unpack climb_rate get_paths t_wind_field
#execs=t_snapshots t_vel_stats m_histogram tracer gl_test gpr_test gpr_test2 gpr_test3 gpr_test4 gpr_test5 ttt_play game_play game_play2 game_play3 game_kf game_frozen game_random interp lp_test random_test t_correl matrix_update var_test t_mr_test

all:
	$(MAKE) executables

executables: $(execs)

depend: $(src)
	$(cxx) $(iflags) -MM $(src) >Makefile.dep

libsrsim.a: $(objs)
	rm -f $@
	ar rs $@ $^

-include Makefile.dep

unpack: unpack.cc common.o
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

get_paths: get_paths.cc common.o
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

t_correl: t_correl.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

climb_rate: climb_rate.cc common.o fileinfo.o glider.o
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

kf_calc: kf_calc.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

t_mr_test: t_mr_test.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

t_snapshots: t_snapshots.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

t_wind_field: t_wind_field.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

t_vel_stats: t_vel_stats.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

m_histogram: m_histogram.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

tracer: tracer.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

interp: interp.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

ttt_play: ttt_play.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

game_play: game_play.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

game_play2: game_play2.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

game_play3: game_play3.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

game_kf: game_kf.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

gpr_la_test: gpr_la_test.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

lapack_test: lapack_test.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

game_frozen: game_frozen.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

game_random: game_random.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

soar: soar.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

gl_test: gl_test.cc glider_test.o libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

gl_conv: gl_conv.cc glider_test.o libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

gpr_test: gpr_test.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lp_lflags)

gpr_test2: gpr_test2.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

gpr_test3: gpr_test3.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

gpr_test4: gpr_test4.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

gpr_test5: gpr_test5.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

lp_test: lp_test.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

random_test: random_test.cc libsrsim.a
	 $(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

matrix_update: matrix_update.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

var_test: var_test.cc libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -c $<

clean:
	rm -f $(execs) $(objs) libsrsim.a

.PHONY: clean all executables depend
