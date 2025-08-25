# Bayesian Network Meta-regression Models for Multivariate Aggregate Responses with Partially Observed or Completely Missing Within-Treatment Sample Covariance Matrices

## Data
- raw data: MMNAC.csv
- processed data for input: MMNAC.txt

## Code
- main program: MNMA_sim_v2.f
- MCMC procedure: gibbsMNMA_sim_v2.f
- helper functions and subroutines:
  - optim1.f: Nelder & Mead simplex algorithm for function minimization
  - hpd.f: computing 100(1-alpha)% HPD and credible intervals
  - utility.f: include all other functions and subroutines used in the program

## Instruction

The code is written in Fortran 90, which requires a Fortran compiler. After loading the compiler, the user must run compile.sh to compile the code. 
