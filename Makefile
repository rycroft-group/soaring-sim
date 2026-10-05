# Load the common configuration file
ifneq ($(wildcard ../config.mk),)
    include ../config.mk
else ifneq ($(wildcard config.mk),)
    include config.mk
else
    $(error No config.mk found. For a stand-alone build, copy a template, e.g. "cp config/config.mk.linux config.mk")
endif

iflags=-Isrc `gsl-config --cflags` $(fftw_iflags)
lflags=`gsl-config --libs` $(fftw_lflags)

# Library source code files
obj_b=turb_fluid.o tf_grid.o tf_grid_mr.o glider.o glider_test.o common.o \
     k_func.o en_spec.o gpr.o mcts.o fileinfo.o soaring_sim.o soaring_sim_io.o wind_correl.o
objs=$(addprefix src/,$(obj_b))

# Main executables
execs=soar unpack climb_rate get_paths

# Additional tests
test_b=en_spectrum gl_conv gl_test kf_calc lapack_test t_correl t_mr_interp \
	   t_mr_test t_vel_stats t_wind_field t_wind_interp w_rms_modes
tests=$(addprefix tests/,$(test_b))

all: $(execs) $(tests)

# Library of all pre-compiled object files from 'src' directory
src/libsrsim.a: $(objs)
	rm -f $@
	ar rs $@ $^

# Main executables
unpack: src/unpack.cc src/common.o
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

get_paths: src/get_paths.cc src/common.o
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

climb_rate: src/climb_rate.cc src/common.o src/fileinfo.o src/glider.o
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags)

soar: src/soar.cc src/libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

# Additional test executables
tests/%: tests/%.cc src/libsrsim.a
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

tests/lapack_test: tests/lapack_test.cc
	$(cxx) $(cflags) $(iflags) -o $@ $^ $(lflags) $(lp_lflags)

%.o: %.cc
	$(cxx) $(cflags) $(iflags) -MMD -MP -c $< -o $@

-include $(objs:.o=.d)

clean:
	rm -f $(execs) $(tests) $(objs) $(objs:.o=.d) src/libsrsim.a

.PHONY: clean all
