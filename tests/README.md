# Soaring-Sim: Additional test programs
This directory contains a number of programs that were used to test and check
different aspects of the soaring simulation code. (See the [top-level
README](../README.md) for a general overview.) A number of these programs were
used to generate data for the main paper and the appendices. The programs are
compiled by typing `make` in the top-level directory. They can then be run from
within this directory. Some programs write their output within this directory,
whereas others write their results to standard output.

## List of programs
- **gl_test** uses the `glider_test` class to numerically integrate a single
  glider over a duration, in a deterministic setting with no wind. The program
  outputs a text file **gt.dat** containing the glider trajectory.

- **gl_conv** uses the `glider_test` class to examine the convergence of the
  glider numerical integration scheme. It first integrates the glider over
  a fixed duration using a very small timestep, and uses that as a reference.
  It then integrates the glider using a range of timesteps *&Delta;t*, and
  outputs the two-norm difference *E* between the solution and the reference.
  Plotting *E* as a function of *&Delta;t* can be used to assess the numerical
  convergence rate. By default the code simulates a glider moving
  deterministically in zero wind, but the effects of a stochastic turbulent
  wind field can be added.

- **kf_calc** tests the computation of the Gaussian process regression kernel.
  It outputs the 2D precomputed table of *K*(*r*,*t*) to the file **KF.dat** in
  a binary format that can be read by the freeware plotting program
  [Gnuplot](https://gnuplot.info). (See the
  [utils-gp](https://github.com/chr1shr/utils-gp)
  documentation for more detail on this format.) To plot the results, use
  the Gnuplot commands
  ```Gnuplot
  set xlabel 'r'
  set ylabel 't'
  set zlabel 'K(r,t)'
  splot 'KF.dat' matrix binary w l
  ```
  The **kf_calc** program also computes a table of a reduced *K*(*r*) for a
  frozen wind field. This is saved to the file **KFr.dat** as plain text and
  can be plotted in Gnuplot using
  ```Gnuplot
  set xlabel 'r'
  set ylabel 'K(r)'
  plot 'KFr.dat' w l
  ```

- **lapack_test** demonstrates how to use the LAPACK library to solve a linear
  system *Ax* = *b* for a 3&times;3 matrix *A*. The program also shows how the
  array holding *A* gets modified during the solve.

- **en_spectrum** initializes a turbulent wind field, and then computes the
  energy spectrum *E*(*k*) as a function of wavenumber *k*. The program can be
  run with the command
  ```Shell
  ./en_spectrum <modes>
  ```
  where `<modes>` is the number of Fourier modes in each direction in the wind
  field. The modes are initialized randomly according to their steady-state
  distributions. An energy spectrum calculation is performed, using a nonlinear
  binning strategy where the *k* values are mapped with a hyperbolic arcsine;
  this results in linear binning for small *k* transitioning to logarithmic
  binning for large *k*. The program time-integrates the modes forward
  following the Ornstein&#8211;Uhlenbeck process described in the paper, and
  performs eight more energy spectrum calculations spaced 100&nbsp;s apart. The
  results, in simulation units, are saved to a file **en_spec_\<modes\>.dat**.
  The first column is the bin index, and the spectrum data is provided in nine
  (*k*,*E*) pairs. If a bin is empty then the pair is written as `NaN NaN`. The
  results confirm that the mode statistics follow Kolmogorov scaling, *E*(*k*)
  &propto; *k*<sup>-5/3</sup> as a function of wave number *k*. In addition,
  the results show that the mode statistics remain invariant as the wind field
  is time-integrated.

- **t_correl** empirically computes the correlation function *K*(*r*,*t*).
  The program requires two arguments and can be run with the command
  ```Shell
  ./t_correl <modes> <max_radius>
  ```
  where `<modes>` is the number of Fourier modes in each direction in the wind
  field and `<max_radius>` is the maximum radius *r* over which to compute
  *K*. Many other parameters are hard-coded into the example, and generally
  match typical parameters used within the paper. The program creates a large
  number of randomly-distributed sample points throughout the turbulent wind
  field domain and evaluates the velocities there to compute *K*. The results
  are outputted to a file **tcor_\<modes\>_r\<radius\>.dat** and can be
  compared to the analytical correlation function.

- **t_mr_test** examines how the wind fields mean-revert over time. A reference
  turbulent wind field is initialized with the random mode distribution in
  statistical steady state. A large number of instances are then simulated,
  starting from this reference state and time-integrating forward by a fixed
  duration. The average wind field at the end of the duration is then stored.
  This is shown to be a close match to deterministically applying exponential
  mean-reversion to all of the modes. Results are all stored within the
  **t_mr_test.odr** directory.

- **t_mr_interp** demonstrates the `turb_fluid_grid_mr` class, which can
  rapidly evaluate the expected wind velocity at a location and time in the
  future, based on mean-reversion of the modes. The class uses a combination of
  trilinear interpolation in space and Hermite interpolation in time as
  described in the paper. For a single location and a range of times into the
  future, the code compares the output of the `turb_fluid_grid_mr` class with a
  direct calculation. Results are stored in the file **t_mr_interp.dat**.

- **t_vel_stats** simulates a turbulent wind field over a duration, and outputs
  statistics about the mean and root-mean-squared (RMS) of the velocity
  components to a file **vstats.dat**. The code computes a long-term average of
  the RMS wind speed. This value can be compared to the imposed RMS wind speed,
  and should match up to stochastic fluctuations.

- **t_wind_field** computes the mode spectrum prefactor *&alpha;* so the RMS
  wind velocity component can match a specific value. The program then outputs
  a sequence of *z* cross-sections of an instance of a turbulent wind field
  with the specified RMS wind velocity component. Results are stored within the
  **twf.odr** output directory.

- **t_wind_interp** compares the grid interpolation functions to an exact
  expansion of the Fourier modes. The general program syntax is
  ```Shell
  ./t_wind_interp <modes> <ax> <ay> <az> <bx> <by> <bz> <steps>
  ```
  where `<modes>` is the number of Fourier modes in each direction in the wind
  field. The program uses a sample line from (`<ax>`,`<ay>`,`<az>`) to
  (`<bx>`,`<by>`,`<bz>`) divided into `<steps>` intervals. The program can
  also be run in a default mode with
  ```Shell
  ./t_wind_interp <modes>
  ```
  in which case it uses a sample line along the *x* axis with
  `<modes>`&times;16 steps. At each point along the sample line, the program
  computes the exact expansion of the wind field in terms of Fourier modes, as
  well as the trilinear, Lanczos-2, and tricubic interpolants. The results are
  outputted to the file **twi.dat**. In addition, if the default mode is used,
  a second file **twi_grid.dat** is outputted that includes the wind velocity
  at the aligned grid points along the *x*-axis.

- **w_rms_modes** calculates how the RMS of a wind velocity component changes
  as a function of the number of modes *N*. As described in the paper, the
  value of the mode spectrum parameter *&alpha;* is fixed so that for *N* = 64
  modes, the RMS of a velocity component is set to 1&nbsp;m/s. This program
  considers a range of *N* from 10 to 160 and computes the RMS of a velocity
  component. Increasing *N* results in a small increase in the RMS of a
  velocity component since there is an extra contribution from the
  short-wavelength modes that were added.
