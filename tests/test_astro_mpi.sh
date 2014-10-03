#!/bin/bash
time -p /usr/bin/mpirun -np 4 ../src/astrochem_mpi -m 1  ./input.ini 
