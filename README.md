# leverage_cycles
This repository contains the code to simulate leverage cycles used to produce my paper "Taming the leverage cycle" (http://arxiv.org/abs/1507.04136).

The code has been tested on Mac OS Yosemite.

Mac OSX Requirements:

- R
- Rcpp (R package for integration with c++)
- R package ggplot2
- C++ compiler

R package installation instructions:

In the folder containing this repository, run the following in the terminal:

1. R CMD INSTALL leverageCycles
2. R CMD check leverageCycles
3. R CMD build leverageCycles

To run the leverage cycle model run test_lev_cycles.R 






