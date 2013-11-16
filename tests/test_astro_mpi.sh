#!/bin/bash
time -p /usr/bin/mpirun -np `cat $OAR_FILE_NODES|wc -l` -machinefile $OAR_FILE_NODES ../src/astrochem  ./input.ini 10
