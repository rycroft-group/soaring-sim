# Soaring-Sim: Library source code and main executables
This directory contains all of the library source code, which is structured
around a number of C++ classes. Four executable programs are also provided.

The main executable is **soar**, which can run a glider soaring simulation, via
a configuration file to specify all of the parameters. The **soar** program can
save a variety of different output files, including the glider trajectories in
a binary format. The **unpack** utility can process these binary files and
output the trajectories as plain text for plotting. Both **soar** and
**unpack** are described in detail in the [top-level README](../README.md).
The other two executables are described here.

## **get_paths**: a utility for outputting Monte Carlo tree search paths
The **soar** simulation has an option for outputting the MCTS paths that
were used for planning. This is enabled with the `path_interval <n>` option,
which causes the paths to be saved after every *n* planning steps.
The paths are stored in a binary file of the form `path_info.<m>` within the
output directory, for the wind field *m*. These binary files can be very
large since they typically contain hundreds of trajectories for each glider.

The **get_paths** utility processes the binary file and extracts one or more
trajectories. The syntax is
```Shell
./get_paths <output_dir> <trial> <glider> <step> <path> <output_file>
```
`<output_dir>` sets the output directory, and `<trial>` sets which wind field
to use. `<glider>` sets the glider number, and `<step>` sets which outputted
planning step to consider.

`<path>` chooses which MCTS path to output. If `<path>` is set to `-`, then all
paths are outputted. `<output_file>` sets where the data is stored, and if `-`
is used then the data is sent to standard output.

Each line in the output file has the format
```
<frame> <sim_time> <path_info> [<path_info> ...]
```
where `<frame>` is an integer output frame number, starting from zero.
`<sim_time>` gives the current time in simulation units, measured from the
start of the simulation.

If `<path>` was set to an integer then there will just be one `<path_info>`
block, but if the `-` option was used then `<path_info>` data for all MCTS
trials will be given. Each `<path_info>` block contains seven entries and has
the form
```
<rx> <ry> <rz> <wx> <wy> <wz> <bank>
```
where the glider's position is (*rx*,*ry*,*rz*) and the glider's prediction of
the wind field is (*wx*,*wy*,*wz*). `<bank>` gives the glider's bank angle
setting, measured as an integer.

## **climb_rate**: a utility for outputting individual climbs
The `summary climb_stats` file described in the [top-level
README](../README.md) provides statistics about the climb rates taken over the
entire set of gliders and wind fields. This file also includes information
about rates of energy gain, taking into account the glider and wind velocity.

The **climb_rate** utility can analyze the minimal binary files
**glider_mb.0**, **glider_mb.1**, ... containing the glider trajectories, and
compute per-glider climb rates. It is intended to be used to understand
detailed aspects of the climb rate distributions. The utility syntax is
```Shell
./climb_rate <input_file> <t_lo> [<t_hi>]
```
The `<input_file>` here is the simulation configuration file, with a **.cfg**
suffix. The utility first reads the configuration file to gather the
conversions from simulation to physical units. `<t_lo>` and `<t_hi>` then set
the time window [*t*<sub>lo</sub>,*t*<sub>hi</sub>] in physical units over
which to compute the climb rates. If `<t_hi>` is omitted then *t*<sub>hi</sub>
is set to the end time of the simulation. The code analyzes the minimal binary
files stored within the output directory (found by replacing the **.cfg**
suffix with the **.odr** suffix).

Output data is stored to a text file **cstats_indiv** within the output
directory. Each line in the file has the form
```
<wind_field_num> <glider_num> <a> <b> <simple_c_rate>
```
The `<wind_field_num>` and `<glider_num>` are integers for the wind field
number and glider number, respectively. Each glider's height *z*(*t*)
is fit using linear regression to *z* = *a* + *bt* and the *a* and *b* values
are shown in the third and fourth entries, respectively; *b* therefore gives a
measure of climb rate. The fifth entry shows an alternative simple climb rate
calculation based on
(*z*(*t*<sub>hi</sub>)-*z*(*t*<sub>lo</sub>))/(*t*<sub>hi</sub>-*t*<sub>lo</sub>).

Note that these computations do not provide information about the rate of
energy gain, since that requires knowledge of the wind field, which is not
included in the trajectory data.

## C++ code structure
A brief overview of the different C++ classes in the library is given below.
The code for a given class is stored within the corresponding files with
**.hh** and/or **.cc** suffixes, unless otherwise noted.

- `glider` – This class represents the state of a single glider, including
  its position, velocity, and bank angle. The class contains routines for
  integrating the glider state via the forward Euler method or improved Euler
  method.

- `glider_model` (in **glider.hh** and **glider.cc**) – This class contains all
  of the glider model parameters. It also calculates several trigonometrical
  tables of bank angles that are used in the differential equations for the
  glider. Only one instance of this class is needed, and can be shared amongst
  many `glider` instances.

- `glider_test` – This class is primarily designed for testing the numerical
  integration of a glider model. The class can be coupled to a turbulent wind
  field, but can also work in a deterministic environment where the wind field is
  assumed to be zero everywhere.

- `turb_fluid` – This class represents a model of three-dimensional
  turbulent wind field made up of a stochastically evolving set of Fourier
  modes. It contains routines for stepping the modes forward in time, as well
  as evaluating the wind at one or several locations.

- `turb_fluid_grid` (in **tf_grid.hh** and **tf_grid.cc**) – This class is
  derived from `turb_fluid`. It contains additional functionality to evaluate
  the wind on a full three-dimensional grid using the FFTW library. It also
  contains routines for rapidly evaluating the wind at any location using
  either trilinear, Lanczos-2, or tricubic interpolation.

- `turb_fluid_grid_mr` (in **tf_grid_mr.hh** and **tf_grid_mr.cc**) – This
  class is derived from `turb_fluid_grid`. It contains additional functionality
  to estimate the wind field into the future assuming that the Fourier modes
  undergo reversion to the mean, following the OU process.

- `kernel_func`, `kernel_rt`, &amp; `kernel_r` (in **k_func.hh** and
  **k_func.cc**) – These classes are used to evaluate the kernel function used
  in the Gaussian process regression (GPR). The kernel function is derived from
  the correlations in the wind field. Since the kernel function needs to be
  called often, the class creates a lookup table of values that are later used
  to perform bilinear interpolation.

- `fileinfo` – This class parses the text configuration files and reads in
  all of the required parameters. It also performs calculations for
  initializing the simulation, such as calculating the timestep.

- `cli_stats` &amp; `mti_stats` (in **stats.hh**) – Small structures for
  computing the mean, standard deviation, minimum, and maximum of a group of
  numbers. `cli_stats` is used for collecting statistics on climb rates, and
  `mti_stats` is used for collecting statistics on the MCTS.

- `wind_correl` – This class is used to compute the correlation between the
  true turbulent wind field and that predicted by Gaussian process regression.
  The class makes use of a small helper structure called `p_cor_coeff` for
  tracking the first and second moments needed to compute the correlation.

- `en_spec_param` (in **en_spec.hh** and **en_spec.cc**) – This class can be
  used to compute the energy spectrum *E*(*k*) as a function of wavenumber *k*
  for the turbulent fluid. The class uses a hyperbolic arcsine
  transformation of the *k* values, which gives roughly linear bins for small
  *k* and roughly logarithmic bins for large *k*. A small helper data structure
  `en_spec_data` is also defined, for storing the energy spectrum data in each
  bin.

- `gpr` – This class can perform Gaussian process regression. It takes
  in a sequence of glider positions and times, along with a corresponding
  measurement of the wind field. Using the kernel function described by
  `kernel_func` and the related classes, the code performs Gaussian process
  regression to predict the wind field at other locations. The class uses the
  Woodbury formula to accelerate the linear algebra computations. It can also
  use the LAPACK library for diagnostic purposes.

- `mcts` – This class performs the Monte Carlo tree search for the glider.
  It builds a tree data structure of possible moves that the glider could make
  along with their favorability for gaining energy. It allows the glider to
  select the action at each step that is likely to lead to the most energy
  gain.

- `soaring_sim` (main code in **soaring_sim.hh** and **soaring_sim.cc**;
  additional output routines in **soaring_sim_io.cc**, and `mcts` connecting
  code in **mcts_sp.cc**). This is the main class for performing the soaring
  simulations. It is derived from the `fileinfo` class, so that it can read in
  all of the required parameters. It makes use of all other classes as
  components, and also contains a number of output routines.

The classes and their main functions are commented in the style of
[Doxygen](https://www.doxygen.nl/), using special comment blocks that begin
with `/**`.
