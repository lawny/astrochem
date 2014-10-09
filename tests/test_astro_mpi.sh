#!/bin/bash
time -p $MPIRUN -np 4 ../src/astrochem_mpi -m 1  ./input.ini 
