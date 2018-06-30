#!/bin/bash 


# OpenFOAM
export PATH=/usr/lib64/mpi/gcc/openmpi/bin/:$PATH
source $HOME/OpenFOAM/OpenFOAM-4.1/etc/bashrc WM_LABEL_SIZE=64

case_dir=$1
solver=$2

blockMesh -case $case_dir &> $case_dir/blockMesh.log

if [ "$solver" = "icoFoam" ]; then
    icoFoam -case $case_dir &> $case_dir/icoFoam.log
elif [ "$solver" = "pimpleDyMFoam" ]; then
    pimpleDyMFoam -case $case_dir &> $case_dir/pimpleDyMFoam.log
elif [ "$solver" = "pisoFoam" ]; then
    pisoFoam -case $case_dir &> $case_dir/pisoFoam.log
else
    echo wrong solver specified! > solver.log
fi

