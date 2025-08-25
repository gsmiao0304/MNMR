#!/bin/bash
module load GCC/13.3.0

# Compile each source file
echo "Compiling..."
gfortran -ffree-form -c MNMA_sim_v2.f || exit 1
gfortran -ffree-form -c gibbsMNMA_sim_v2.f || exit 1
gfortran -ffixed-form -std=legacy -c optim1.f || exit 1
gfortran -ffixed-form -std=legacy -c utility.f || exit 1
gfortran -ffixed-form -std=legacy -c hpd.f || exit 1

# Link
gfortran MNMA_sim_v2.o gibbsMNMA_sim_v2.o optim1.o utility.o hpd.o -o mnma_model || exit 1
echo "Compilation successful!"