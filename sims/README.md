# Soaring-Sim: Sample configuration files
This directory contains several sample configuration files for the soaring
simulation. (See the [top-level README](../README.md) for a general overview.)
Each configuration file has a **.cfg** suffix, and can be passed to the main
simulation executable **soar** with the syntax
```Shell
OMP_NUM_THREADS=<n> ./soar <config_file>
```
where `<n>` is the number of OpenMP threads to use. The results will be stored
in a directory matching the name of the config file, but with the suffix
changed to **.odr**.

- **lfd64small.cfg** – This file is provided as a small initial example, and
  simulates 4 turbulent wind fields, each with 8 gliders. The simulation uses
  *&tau;<sub>c</sub> = 64&nbsp;s to set the temporal correlations of the
  stochastically evolving wind Fourier modes. the full information glider model,
  where the Monte Carlo tree search planning assumes that the glider can see
  the complete wind field at the current time. The simulation should take
  around 10 min to run, and creates much less data than most trials performed
  for the paper.

- **gd48m20.cfg** – This file simulates partial information gliders in wind
  fields with *&tau;<sub>c</sub> = 32&nbsp;s. The gliders build a model of the
  wind using Gaussian process regression (GPR) with a memory duration of
  20&nbsp;s. The file simulates 32 wind fields, each with 24 gliders; this is
  a standard 

- **gd48m20d6.cfg** – This is a variant of the **gd24m20d6.cfg** where the MCTS
  planning depth is changed from 12 to 6. As described in the paper this
  results in a substantial loss of performance.

- **id64N16.cfg** – This configuration file simulates large sample of 2048
  wind fields, each with 24 gliders. It uses *&tau;<sub>c</sub> = 64&nbsp;s,
  and the full information glider model. The wind field uses *N* = 16 Fourier
  modes in each direction. This file was used to analyze climb rate
  distributions, and search for outliers. (Note that running this file can take
  several days of computation time.)

- **rd24.cfg** – This file simulates gliders moving randomly in wind fields
  with *&tau;<sub>c</sub> = 24&nbsp;s. The random glider climb rates were used
  in the paper as a lower baseline on glider performance.

- **zdf.cfg** – This file simulates gliders flying straight in a frozen wind
  field. Gliders flying straight were used as a point of comparison in the
  paper. 
