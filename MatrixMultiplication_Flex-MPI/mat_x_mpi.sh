#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

DIR="$(dirname -- $0)"

port1=6666
port2=6667
nodefile=${DIR}/mat_x_mpi_nodefile
executable=${DIR}/mat_x_mpi
controllernode=slurm-node-1
rankfile=${DIR}/mat_x_mpi_rankfile

#executin m_wacomm1 program
mpiexec -genvall -np $1 -f $rankfile  ${executable} -cfile ${nodefile} -policy-malleability-triggered -lbpolicy-static -ni 2 -ports $port1 $port2 -controller ${controllernode} -IOaction 2 -alloc:0
